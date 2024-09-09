/** Copyright Stewart Allen <sa@grid.space> -- All Rights Reserved */

"use strict";

// toast/burnt-edge material and helpers

// modified EdgesGeometry to preserve edge/face relationships
function createEdgeData(obj, thresholdAngle = 20) {

    const { geometry, position } = obj;

    const { MathUtils, Triangle, Vector3, BufferAttribute, ShaderMaterial } = THREE;
    const _v0 = /*@__PURE__*/ new Vector3();
    const _v1 = /*@__PURE__*/ new Vector3();
    const _normal = /*@__PURE__*/ new Vector3();
    const _triangle = /*@__PURE__*/ new Triangle();

    const precisionPoints = 4;
    const precision = Math.pow(10, precisionPoints);
    const thresholdDot = Math.cos(MathUtils.DEG2RAD * thresholdAngle);
    let positionAttr = geometry.getAttribute('position');

    const decim = base.verticesToPoints(positionAttr.array, {
        threshold: 1,
        precision: 0.05,
        maxpass: 10
    }).map(p => p.toArray()).flat();

    console.log({ decim, array: positionAttr.array });
    positionAttr = new BufferAttribute(decim.toFloat32(), 3);
    geometry.setAttribute('position', positionAttr);

    const posArray = positionAttr.array;

    const indexCount = posArray.length;
    const numFaces = indexCount / 3;
    const indexArr = [0, 0, 0];
    const vertKeys = ['a', 'b', 'c'];
    const hashes = new Array(3);
    const edgeData = {};
    const vertices = [];
    const faceRecords = new Array(numFaces); // [i1, i2] array per face
    const pointToEdge = {};

    function update(rec, val) {
        let arr = rec.idx;
        if (arr.indexOf(val) >= 0) {
            return;
        }
        let io = arr.indexOf(0);
        if (io >= 0) {
            arr[io] = val;
            rec.match++;
        } else {
            console.log({ arr, io, val });
            throw "array full";
        }
    }

    function updatePoints(pointKey, lineIndex) {
        let edges = pointToEdge[pointKey];
        if (!edges) {
            edges = (pointToEdge[pointKey] = []);
        }
        edges.addOnce(lineIndex);
    }

    for (let i = 0, faceId = 0; i < indexCount; i += 3, faceId++) {
        // points from face
        indexArr[0] = i;
        indexArr[1] = i + 1;
        indexArr[2] = i + 2;

        const { a, b, c } = _triangle;
        a.fromBufferAttribute(positionAttr, indexArr[0]);
        b.fromBufferAttribute(positionAttr, indexArr[1]);
        c.fromBufferAttribute(positionAttr, indexArr[2]);
        _triangle.getNormal(_normal);

        // create point hashes for the edge from the vertices
        hashes[0] = `${Math.round(a.x * precision)},${Math.round(a.y * precision)},${Math.round(a.z * precision)}`;
        hashes[1] = `${Math.round(b.x * precision)},${Math.round(b.y * precision)},${Math.round(b.z * precision)}`;
        hashes[2] = `${Math.round(c.x * precision)},${Math.round(c.y * precision)},${Math.round(c.z * precision)}`;

        // create face record
        let rec = faceRecords[faceId] = { points:hashes.slice(), idx:[0,0,0,0], match:0 };

        // skip degenerate triangles
        if (hashes[0] === hashes[1] || hashes[1] === hashes[2] || hashes[2] === hashes[0]) {
            continue;
        }

        // iterate over every edge
        for (let j = 0; j < 3; j++) {
            // get the first and next vertex making up the edge
            const jNext = (j + 1) % 3;
            const vecHash0 = hashes[j];
            const vecHash1 = hashes[jNext];
            const v0 = _triangle[vertKeys[j]];
            const v1 = _triangle[vertKeys[jNext]];
            const hash = `${vecHash0}_${vecHash1}`;
            const reverseHash = `${vecHash1}_${vecHash0}`;

            const adjacent = edgeData[reverseHash];
            if (adjacent) {
                // if we found a sibling edge add it into the vertex array if
                // it meets the angle threshold and delete the edge from the map.
                if (_normal.dot(edgeData[reverseHash].normal) <= thresholdDot) {
                    let lineIndex = vertices.length / 4 + 1;
                    vertices.push(v0.x, v0.y, v0.z, 0);
                    vertices.push(v1.x, v1.y, v1.z, 0);
                    // add line index to face and adjoining face
                    update(rec, lineIndex);
                    update(faceRecords[adjacent.faceId], lineIndex);
                    // update points on line records
                    updatePoints(vecHash0, lineIndex);
                    updatePoints(vecHash1, lineIndex);
                }
                edgeData[reverseHash] = null;
            } else if (!(hash in edgeData)) {
                // if we've already got an edge here then skip adding a new one
                edgeData[hash] = {
                    index0: indexArr[j],
                    index1: indexArr[jNext],
                    normal: _normal.clone(),
                    faceId
                };
            }
        }
    }

    // iterate over all remaining, unmatched edges and add them to the vertex array
    if (false)
    for ( const key in edgeData ) {
        if ( edgeData[ key ] ) {
            const { index0, index1 } = edgeData[ key ];
            _v0.fromBufferAttribute( positionAttr, index0 );
            _v1.fromBufferAttribute( positionAttr, index1 );
            vertices.push( _v0.x, _v0.y, _v0.z, 0 );
            vertices.push( _v1.x, _v1.y, _v1.z, 0 );
            console.log('missed', key);
        }
    }

    // for faces with no edges matches to lines, check if any of their
    // points is matched with a line and add those
    for (let rec of faceRecords.filter(r => r.match === 0)) {
        let faces = rec.points.map(pk => pointToEdge[pk]);
        // console.log(faces);
    }

    // Create a BufferAttribute for the line indices
    const numVertices = vertices.length / 4;  // Assuming each line point has 4 components (x, y, z, w)
    const textureDimension = Math.ceil(Math.sqrt(numVertices));  // Get a near-square dimension
    const totalEdges = numVertices / 2;
    // Adjust the vertices array to fit into a 2D texture
    const paddedLength = textureDimension * textureDimension * 4;  // 4 components per vertex (RGBA)
    const paddedVertices = new Float32Array(paddedLength);
    // Copy the original vertices into the padded array
    paddedVertices.set(vertices.toFloat32());

    const edgeLines = new THREE.DataTexture(paddedVertices, textureDimension, textureDimension, THREE.RGBAFormat, THREE.FloatType);
    edgeLines.needsUpdate = true;

    const material = new ShaderMaterial({
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
        uniforms: {
            burnRadius: { value: 1.0 },
            edgeLines: { value: edgeLines },
            totalEdges: { value: totalEdges },
            textureDimension: { value: textureDimension },
        }
    });

    const mgeo = geometry.clone();
    const mesh = new THREE.Mesh(mgeo, material);

    console.log({ numFaces, faceRecords, pointToEdge, edges: vertices.group(4) });

    return { vertices, material, mesh };
}

const vertexShader = `
precision highp float;
precision highp int;
#include <common>
#include <logdepthbuf_pars_vertex>

uniform sampler2D edgeLines;     // Texture containing edge endpoints
varying vec3 vPosition;

void main() {
    vPosition = (modelMatrix * vec4(position, 1.0)).xyz;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    #include <logdepthbuf_vertex>
}
`;

const fragmentShader = `
precision highp float;
precision highp int;
#include <common>
#include <logdepthbuf_pars_fragment>

uniform sampler2D edgeLines;      // Texture containing edge line endpoints
uniform mat4 modelMatrix;         // Pass the modelMatrix to the fragment shader as a uniform
uniform int totalEdges;           // Total unpadded line count
uniform float burnRadius;         // Proximity radius for burn effect
uniform float textureDimension;   // The width and height of the square texture
varying vec3 vPosition;           // Vertex position passed from the vertex shader

// Fetch a line segment from the edgeLines texture and transform to world space
vec3 getEdgeLine(int index) {
    int texWidth = int(textureDimension);                                         // Texture width
    int row = index / texWidth;                                                   // Row in the texture
    int col = index % texWidth;                                                   // Column in the texture
    vec2 uv = vec2(float(col) / textureDimension, float(row) / textureDimension); // Normalize the texture coordinates (UV)
    vec3 linePoint = texture2D(edgeLines, uv).rgb;                                // Fetch the line point from the 2D texture
    return (modelMatrix * vec4(linePoint, 1.0)).xyz;                              // Transform the line position to world space using modelMatrix
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

    // Loop through all line segments in the edgeLines texture
    for (int i = 0; i < totalEdges; i++) {
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
`;

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
        create: createEdgeData
    });
});
