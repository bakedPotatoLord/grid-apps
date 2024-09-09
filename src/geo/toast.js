/** Copyright Stewart Allen <sa@grid.space> -- All Rights Reserved */

"use strict";

// toast/burnt-edge material and helpers

// modified EdgesGeometry to preserve edge/face relationships
function createedgeHash(obj, thresholdAngle = 20) {

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
        pointToEdge = {};

    function updatePoint(pointKey, lineIndex) {
        let edges = pointToEdge[pointKey];
        if (!edges) {
            edges = (pointToEdge[pointKey] = []);
        }
        edges.addOnce(lineIndex);
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

        // create point hashes for the edge
        hashes[0] = `${Math.round(a.x * precision)},${Math.round(a.y * precision)},${Math.round(a.z * precision)}`;
        hashes[1] = `${Math.round(b.x * precision)},${Math.round(b.y * precision)},${Math.round(b.z * precision)}`;
        hashes[2] = `${Math.round(c.x * precision)},${Math.round(c.y * precision)},${Math.round(c.z * precision)}`;

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
                    rec.points.remove(vecHash0);
                    rec.points.remove(vecHash1);
                    // update adjacent face
                    let adj = faces[adjacent.faceId];
                    adj.edges.push(lineIndex);
                    adj.points.remove(adjacent.vecHash0);
                    adj.points.remove(adjacent.vecHash1);
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

    // iterate over all remaining, unmatched edges and add them to the vertex array
    if (false)
    for ( let key in edgeHash ) {
        if ( edgeHash[ key ] ) {
            let { index0, index1 } = edgeHash[ key ];
            _v0.fromBufferAttribute( positionAttr, index0 );
            _v1.fromBufferAttribute( positionAttr, index1 );
            edgeData.push( _v0.x, _v0.y, _v0.z, 0 );
            edgeData.push( _v1.x, _v1.y, _v1.z, 0 );
            console.log('missed', key);
        }
    }

    // for faces with no edges matches to lines, check if any of their
    // points is matched with a line and add those
    for (let rec of faces.filter(r => r.points.length)) {
        let edges = rec.points.map(pk => pointToEdge[pk]).flat();
        rec.edges.push(...edges);
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

    // TODO: create THREE.DataTexture from faceEdgeData

    // Create a BufferAttribute for the line indices
    let edgePoints = edgeData.length / 4;  // Assuming each line point has 4 components (x, y, z, w)
    let edgeDims = Math.ceil(Math.sqrt(edgePoints));  // Get a near-square dimension
    let edgeCount = edgePoints / 2;
    // Adjust the edgeData array to fit into a 2D texture
    let edgesPadLen = edgeDims * edgeDims * 4;  // 4 components per vertex (RGBA)
    let edgeTextData = new Float32Array(edgesPadLen);
    edgeTextData.set(edgeData.toFloat32());

    let edgeLines = new THREE.DataTexture(edgeTextData, edgeDims, edgeDims, THREE.RGBAFormat, THREE.FloatType);
    edgeLines.needsUpdate = true;

    let material = new ShaderMaterial({
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
        uniforms: {
            burnRadius: { value: 1.0 },
            edgeLines: { value: edgeLines },
            edgeCount: { value: edgeCount },
            edgeDims: { value: edgeDims },
            // TODO: pass faceDims and faceIndices
            // TODO: no longer need edgeCount (in shader)
        }
    });

    let mgeo = new THREE.BufferGeometry().setAttribute('position', positionAttr);
    let mesh = new THREE.Mesh(mgeo, material);

    console.log({
        faces,
        faceEdgeData,
        // edges: edgeData.group(4)
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
// TODO: add faceDims and faceIndices

uniform sampler2D edgeLines;      // Texture containing edge line endpoints
uniform mat4 modelMatrix;         // Pass the modelMatrix to the fragment shader as a uniform
uniform int edgeCount;            // Total unpadded line count
uniform float burnRadius;         // Proximity radius for burn effect
uniform float edgeDims;           // The width and height of the square texture
varying vec3 vPosition;           // Vertex position passed from the vertex shader
varying float vFaceIndex;         // Face index passed from the vertex shader

// Fetch a line segment from the edgeLines texture and transform to world space
vec3 getEdgeLine(int index) {
    int texWidth = int(edgeDims);                                             // Texture width
    int row = index / texWidth;                                               // Row in the texture
    int col = index % texWidth;                                               // Column in the texture
    vec2 uv = vec2(float(col) / edgeDims, float(row) / edgeDims); // Normalize the texture coordinates (UV)
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

void main() {
    vec4 color = vec4(1.0, 1.0, 1.0, 1.0);  // Base color (white)
    float burnFactor = 0.0;                 // Accumulate burn factor here

    // TODO: replace with lookup in faceIndices
    // TODO: jump to offset record and read indices until 0
    // TODO: calling getEdgeLine() as below
    // Loop through all line segments in the edgeLines texture
    for (int i = 0; i < edgeCount; i++) {
        vec3 lineStart = getEdgeLine(i * 2);                                          // Fetch the start point of the line segment
        vec3 lineEnd = getEdgeLine(i * 2 + 1);                                        // Fetch the end point of the line segment
        float distToLine = distanceToLine(vPosition, lineStart, lineEnd);             // Calculate the proximity (distance) to the current line segment
        burnFactor = max(burnFactor, clamp(1.0 - distToLine / burnRadius, 0.0, 1.0)); // Linear gradient effect based on distance to line
    }

    #include <logdepthbuf_fragment>

    // Apply darkening effect based on burn factor
    vec4 burntColor = vec4(color.rgb * (1.0 - burnFactor * 0.8), color.a);  // Darken the color near edges
    gl_FragColor = burntColor;
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
