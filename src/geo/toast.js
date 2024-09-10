/** Copyright Stewart Allen <sa@grid.space> -- All Rights Reserved */

"use strict";

// toast/burnt-edge material and helpers

// modified EdgesGeometry to preserve edge/face relationships
function createedgeHash(obj, thresholdAngle = 20, burnRadius = 1) {

    let { geometry } = obj,
        { MathUtils, Triangle, Vector3, BufferAttribute, ShaderMaterial } = THREE,
        _v0 = new Vector3(),
        _v1 = new Vector3(),
        _normal = new Vector3(),
        _triangle = new Triangle(),
        precisionPoints = 4,
        precision = Math.pow(10, precisionPoints),
        thresholdDot = Math.cos(MathUtils.DEG2RAD * thresholdAngle),
        positionAttr = geometry.getAttribute('position');

    if (true) {
        // optional decimation
        let decim = base.verticesToPoints(positionAttr.array, {
            threshold: 1,
            precision: 0.05,
            maxpass: 10
        }).map(p => p.toArray()).flat();
        console.log({ decimation_ratio: decim.length / positionAttr.array.length });
        positionAttr = new BufferAttribute(decim.toFloat32(), 3);
    }

    let posArray = positionAttr.array,
        numVerts = posArray.length / 3,
        numFaces = numVerts / 3,
        vertIndx = [0, 0, 0],
        vertKeys = ['a', 'b', 'c'],
        faces = new Array(numFaces),
        hashes = new Array(3),
        edgeHash = {},
        edgeData = [],
        pointRecs = {};

    function pointLineDistance(point, start, end) {
        const lineDir = new THREE.Vector3().subVectors(end, start);
        const pointToStart = new THREE.Vector3().subVectors(point, start);
        const t = THREE.MathUtils.clamp(pointToStart.dot(lineDir) / lineDir.lengthSq(), 0, 1);
        const projection = new THREE.Vector3().copy(lineDir).multiplyScalar(t).add(start);
        return point.distanceTo(projection);
    }

    function createHash(vec, pos) {
        let key = hashes[pos] = `${Math.round(vec.x * precision)},${Math.round(vec.y * precision)},${Math.round(vec.z * precision)}`;
        let rec = pointRecs[key];
        if (!rec) {
            pointRecs[key] = {
                point: vec.clone(),
                edges: [],
            };
        }
    }

    function updatePoint(pointKey, lineIndex) {
        pointRecs[pointKey].edges.addOnce(lineIndex);
    }

    for (let i = 0, faceId = 0; i < numVerts; i += 3, faceId++) {
        // vertex indices into position array
        vertIndx[0] = i;
        vertIndx[1] = i + 1;
        vertIndx[2] = i + 2;

        // reuse triangle vector3s
        let { a, b, c } = _triangle;
        a.fromBufferAttribute(positionAttr, vertIndx[0]);
        b.fromBufferAttribute(positionAttr, vertIndx[1]);
        c.fromBufferAttribute(positionAttr, vertIndx[2]);

        // compute normal for connected face comparison
        _triangle.getNormal(_normal);

        createHash(a, 0);
        createHash(b, 1);
        createHash(c, 2);

        // create face record
        let rec = faces[faceId] = {
            edges: [],
            points: hashes.slice()
        };

        // skip degenerate triangles
        if (hashes[0] === hashes[1] || hashes[1] === hashes[2] || hashes[2] === hashes[0]) {
            continue;
        }

        // iterate over every edge
        for (let j = 0; j < 3; j++) {
            // get the first and next vertex making up the edge
            let jNext = (j + 1) % 3;
            let vecHash0 = hashes[j];
            let vecHash1 = hashes[jNext];
            let v0 = _triangle[vertKeys[j]];
            let v1 = _triangle[vertKeys[jNext]];

            // adjacent edges should have reversed point order (if manifold)
            let hash = `${vecHash0}_${vecHash1}`;
            let reverseHash = `${vecHash1}_${vecHash0}`;

            let adjacent = edgeHash[reverseHash];
            if (adjacent) {
                // if we found a sibling edge add it into the vertex array if
                // it meets the angle threshold and delete the edge from the map
                if (_normal.dot(edgeHash[reverseHash].normal) <= thresholdDot) {
                    let lineIndex = edgeData.length / 4 + 1;
                    edgeData.push(v0.x, v0.y, v0.z, 0);
                    edgeData.push(v1.x, v1.y, v1.z, 0);
                    // add line index to face and adjoining face
                    rec.edges.push(lineIndex);
                    // update adjacent face
                    let adj = faces[adjacent.faceId];
                    adj.edges.push(lineIndex);
                    // update points on line records
                    updatePoint(vecHash0, lineIndex);
                    updatePoint(vecHash1, lineIndex);
                }
                edgeHash[reverseHash] = null;
            } else if (!(hash in edgeHash)) {
                // if we've already got an edge here then skip adding a new one
                edgeHash[hash] = {
                    index0: vertIndx[j],
                    index1: vertIndx[jNext],
                    normal: _normal.clone(),
                    faceId,
                    vecHash0,
                    vecHash1
                };
            }
        }
    }

    // for points with no associated edges, find edge lines within burnRadius
    for (let [ key, val ] of Object.entries(pointRecs)) {
        for (let i=0; i<edgeData.length; i += 8) {
            _v0.set(edgeData[i + 0], edgeData[i + 1], edgeData[i + 2]);
            _v1.set(edgeData[i + 4], edgeData[i + 5], edgeData[i + 6]);
            let dist = pointLineDistance(val.point, _v0, _v1);
            if (dist <= burnRadius) {
                console.log('ADD', { dist, p:val.point, v0:_v0, v1:_v1} );
                val.edges.addOnce(i/4 + 1);
                // TODO if a point is near an edge, check the min distance
                // between the two lines connected to this point and all edges
            }
        }
    }

    // for faces with no edges matches to lines, check if any of their
    // points is matched with a line and add those]
    let maxEdges = 0;
    let minEdges = Infinity;
    for (let rec of faces) {
        let edges = rec.points.map(pk => pointRecs[pk].edges).flat();
        // console.log(edges);
        // TODO fix ... why do we need at least one value before 0??
        rec.edges = [...rec.edges, ...edges, -1].uniq();
        maxEdges = Math.max(maxEdges, rec.edges.length);
        minEdges = Math.min(minEdges, rec.edges.length);
    }

    // encode face to edge index array buffer
    let faceEdgesOff = [];
    let faceEdgesIdx = [];
    for (let rec of faces) {
        faceEdgesOff.push(faceEdgesIdx.length + numFaces);
        faceEdgesIdx.push(...rec.edges, 0);
    }
    let faceEdgeData = new Float32Array(faceEdgesOff.length + faceEdgesIdx.length);
    faceEdgeData.set(faceEdgesOff);
    faceEdgeData.set(faceEdgesIdx, numFaces);
    // Calculate the size of the texture to store the face-to-edge data
    let faceDims = Math.ceil(Math.sqrt(faceEdgeData.length));
    // Create a padded array to match the texture size
    let faceTextData = new Float32Array(faceDims * faceDims);  // one float per entry
    faceTextData.set(faceEdgeData);  // Copy the faceEdgeData into the padded array
    // Create the DataTexture using the faceTextData
    let faceIndices = new THREE.DataTexture(faceTextData, faceDims, faceDims, THREE.RedFormat, THREE.FloatType);
    faceIndices.needsUpdate = true;

    // Create a BufferAttribute for the line indices
    let edgeDims = Math.ceil(Math.sqrt(edgeData.length / 4));  // Get a near-square dimension
    // Adjust the edgeData array to fit into a 2D texture
    let edgesPadLen = edgeDims * edgeDims * 4;  // 4 components per vertex (RGBA)
    let edgeTextData = new Float32Array(edgesPadLen);
    edgeTextData.set(edgeData.toFloat32());
    // Create the DataTexture using the edgeTextData
    let edgeLines = new THREE.DataTexture(edgeTextData, edgeDims, edgeDims, THREE.RGBAFormat, THREE.FloatType);
    edgeLines.needsUpdate = true;

    let material = new ShaderMaterial({
        transparent: true,
        blending: THREE.NormalBlending,
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
        uniforms: {
            burnRadius: { value: burnRadius },
            edgeLines: { value: edgeLines },
            edgeDims: { value: edgeDims },
            faceIndices: { value: faceIndices },
            faceDims: { value: faceDims },
            maxEdges: { value: maxEdges }
        }
    });

    let mgeo = new THREE.BufferGeometry().setAttribute('position', positionAttr);
    let mesh = new THREE.Mesh(mgeo, material);

    console.log({
        minEdges,
        maxEdges,
        pointRecs,
        faces,
        faceDims,
        faceTextData,
        edges: edgeData.group(4),
        edgeDims,
        edgeTextData,
    });

    return { edgeData, material, mesh };
}

let vertexShader = `
precision highp float;
precision highp int;
#include <common>
#include <logdepthbuf_pars_vertex>

uniform sampler2D edgeLines;
varying vec3 vPosition;
varying float vFaceIndex;

void main() {
    vFaceIndex = floor(float(gl_VertexID) / 3.0);
    vPosition = (modelMatrix * vec4(position, 1.0)).xyz;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    #include <logdepthbuf_vertex>
}
`.trim();

let fragmentShader = `
precision highp float;
precision highp int;
#include <common>
#include <logdepthbuf_pars_fragment>
uniform sampler2D faceIndices;    // The DataTexture containing face-to-edge indices
uniform float faceDims;           // The width and height of the face index texture

uniform sampler2D edgeLines;      // Texture containing edge line endpoints
uniform mat4 modelMatrix;         // Pass the modelMatrix to the fragment shader as a uniform
uniform float burnRadius;         // Proximity radius for burn effect
uniform float edgeDims;           // The width and height of the square texture
uniform int maxEdges;             // max edge count from geometry
varying vec3 vPosition;           // Vertex position passed from the vertex shader
varying float vFaceIndex;         // Face index passed from the vertex shader

// Fetch a line segment from the edgeLines texture and transform to world space
vec3 getEdgeLine(int index) {
    int width = int(edgeDims);                                             // Texture width
    int row = index / width;                                               // Row in the texture
    int col = index % width;                                               // Column in the texture
    vec2 uv = vec2(float(col) / edgeDims, float(row) / edgeDims);             // Normalize the texture coordinates (UV)
    vec3 linePoint = texture2D(edgeLines, uv).rgb;                            // Fetch the line point from the 2D texture
    return (modelMatrix * vec4(linePoint, 1.0)).xyz;                          // Transform the line position to world space using modelMatrix
}

// Calculate the distance from the fragment to the closest point on a line segment
float distanceToLine(vec3 point, vec3 start, vec3 end) {
    vec3 lineDir = normalize(end - start);                                     // Direction of the line
    vec3 v = point - start;                                                    // Vector from line start to the point
    float d = dot(v, lineDir);                                                 // Project the point onto the line
    vec3 closestPoint = start + clamp(d, 0.0, length(end - start)) * lineDir;  // Find the closest point on the line
    return length(closestPoint - point);                                       // Return the distance from the point to the closest point on the line
}

int textureFetch(sampler2D textdata, int index, float texDims) {
    // Calculate row and column based on the index
    int texWidth = int(texDims);
    int row = index / texWidth;
    int col = index % texWidth;

    // Compute the UV coordinates for this row and column
    vec2 uv = vec2(float(col) / texDims, float(row) / texDims);

    // Clamp UV to ensure it's within the valid range (0 to 1)
    uv = clamp(uv, vec2(0.0), vec2(1.0));

    // Fetch the value from the texture and return it as an integer
    return int(texture2D(textdata, uv).r);
}

void main() {
    vec4 color = vec4(1.0, 1.0, 1.0, 0.8);  // Base color (white)
    float burnFactor = 0.0;                 // Accumulate burn factor here

    // Fetch the offset in the faceEdgeData (offset to the start of edge indices for this face)
    int faceOffset = textureFetch(faceIndices, int(vFaceIndex), faceDims);

    // Iterate through the edge indices for this face until we encounter a 0
    // Assuming a reasonable maximum number of edges per face
    for (int i = 0; i < maxEdges; i++) {
        int edgeIndex = textureFetch(faceIndices, i + faceOffset, faceDims);
        if (edgeIndex == 0) break;

        // Call getEdgeLine to fetch and process each edge line segment
        vec3 lineStart = getEdgeLine(edgeIndex - 1);
        vec3 lineEnd = getEdgeLine(edgeIndex);
        float distToLine = distanceToLine(vPosition, lineStart, lineEnd);
        burnFactor = max(burnFactor, clamp(1.0 - distToLine / burnRadius, 0.0, 1.0));
    }

    #include <logdepthbuf_fragment>

    // Apply darkening effect based on burn factor
    if (burnFactor > 0.0) {
        vec4 burntColor = vec4(color.rgb * (1.0 - burnFactor * 0.8), 1.0);
        gl_FragColor = burntColor;
    } else {
        gl_FragColor = color;
    }
}
`.trim();

// function addToast(material) {
//     material.onBeforeCompile = (shader) => {
//         shader.vertexShader = shader.vertexShader.replace(
//             `#include <worldpos_vertex>`,
//             `
//             #include <worldpos_vertex>
//             vWorldPosition = vec3(transformed);
//             `
//         );

//         shader.vertexShader = `
//             varying vec3 vWorldPosition;
//         ` + shader.vertexShader;

//         shader.fragmentShader = `
//             varying vec3 vWorldPosition;
//         ` + shader.fragmentShader;

//         shader.fragmentShader = shader.fragmentShader.replace(
//             `#include <dithering_fragment>`,
//             `
//             #include <dithering_fragment>
//             if (vWorldPosition.z < 0.0) {
//                 gl_FragColor.rgb += vec3(0.5, 0.0, 0.0); // Add red tint
//             }
//             `
//         );
//     };
//     return material;
// }

// dep: add.three
// dep: geo.points
// dep: moto.license
gapp.register("geo.toast", [], (root, exports) => {
    exports({
        create: createedgeHash
    });
});
