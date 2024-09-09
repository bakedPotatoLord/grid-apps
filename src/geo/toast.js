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
    const positionAttr = geometry.getAttribute('position');
    const indexCount = positionAttr.count;
    const numFaces = indexCount / 3;
    const indexArr = [0, 0, 0];
    const vertKeys = ['a', 'b', 'c'];
    const hashes = new Array(3);
    const edgeData = {};
    const vertices = [];
    const faceRecords = new Array(numFaces); // [i1, i2] array per face
    const pointToEdge = {};

    function updatePoints(pointKey, lineIndex) {
        let edges = pointToEdge[pointKey];
        if (!edges) {
            edges = (pointToEdge[pointKey] = []);
        }
        edges.addOnce(lineIndex);
    }

    function distToLine(A, B, P) {
        const AB = new THREE.Vector3().subVectors(B, A);   // Vector from A to B
        const AP = new THREE.Vector3().subVectors(P, A);   // Vector from A to P
        const lengthAB = AB.length();  // Length of AB
        const AB_normalized = AB.clone().normalize(); // Normalize AB to unit length
        // Project AP onto AB, this gives us the point on the line closest to P
        const projection = AP.dot(AB_normalized);
        // Clamp the projection to the bounds of the segment
        const t = Math.max(0, Math.min(lengthAB, projection));
        // Find the closest point on the segment to P
        const closestPoint = AB_normalized.multiplyScalar(t).add(A);
        // Return the distance between the closest point and P
        return closestPoint.distanceTo(P);
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
        let rec = faceRecords[faceId] = {
            points:hashes.slice(),
            on:0,
            dist:[
                distToLine(b,c,a),
                distToLine(c,a,b),
                distToLine(a,b,c),
            ],
            len:[
                a.distanceTo(b),
                b.distanceTo(c),
                c.distanceTo(a),
            ]
        };

        // skip degenerate triangles
        if (hashes[0] === hashes[1] || hashes[1] === hashes[2] || hashes[2] === hashes[0]) {
            continue;
        }

        // iterate over every edge
        for (let j = 0; j < 3; j++) {
            // get the first and next vertex making up the edge
            const side = (j + 2) % 3;
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
                    // update points on line records
                    updatePoints(vecHash0, lineIndex);
                    updatePoints(vecHash1, lineIndex);
                    // update "on" flag
                    let adj = faceRecords[adjacent.faceId];
                    rec.on |= (1 << side);
                    adj.on |= (1 << adjacent.side);
                }
                edgeData[reverseHash] = null;
            } else if (!(hash in edgeData)) {
                // if we've already got an edge here then skip adding a new one
                edgeData[hash] = {
                    index0: indexArr[j],
                    index1: indexArr[jNext],
                    normal: _normal.clone(),
                    faceId,
                    side
                };
            }
        }
    }

    const material = new ShaderMaterial({
        transparent: true,
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
        uniforms: {
            burnRadius: { value: 1 },
        }
    });

    const color1 = new BufferAttribute(new Float32Array(numFaces * 3 * 2), 2)
    const array1 = color1.array.fill(1);
    const color2 = new BufferAttribute(new Float32Array(numFaces * 3 * 2), 2)
    const array2 = color2.array.fill(1);
    const color3 = new BufferAttribute(new Float32Array(numFaces * 3 * 2), 2)
    const array3 = color3.array.fill(1);

    for (let i = 0, fid = 0; i < array1.length; i += 6) {
        const rec = faceRecords[fid++];
        const [ a, b, c ] = rec.dist;
        const [ al, bl, cl ] = rec.len;

        const onA = rec.on & 1;
        const onB = rec.on & 2;
        const onC = rec.on & 4;

        console.log(fid, rec.on, rec.len);

        if (onA) {
            array1[i + 0] = a;
            array1[i + 1] = a;
            array1[i + 2] = 0;
            array1[i + 3] = 0;
            array1[i + 4] = 0;
            array1[i + 5] = 0;
        }

        if (onB) {
            array2[i + 0] = 0;
            array2[i + 1] = 0;
            array2[i + 2] = b;
            array2[i + 3] = b;
            array2[i + 4] = 0;
            array2[i + 5] = 0;
        }

        if (onC) {
            array3[i + 0] = 0;
            array3[i + 1] = 0;
            array3[i + 2] = 0;
            array3[i + 3] = 0;
            array3[i + 4] = c;
            array3[i + 5] = c;
        }
    }

    const mgeo = geometry.clone();
    mgeo.setAttribute('color1', color1);
    mgeo.setAttribute('color2', color2);
    mgeo.setAttribute('color3', color3);
    const mesh = new THREE.Mesh(mgeo, material);

    console.log({ numFaces, faceRecords, pointToEdge, edges: vertices.group(4) });

    return { vertices, material, mesh };
}

const vertexShader = `
precision highp float;
precision highp int;

#include <common>
#include <logdepthbuf_pars_vertex>

varying vec3 vPosition;

attribute vec2 color1;
attribute vec2 color2;
attribute vec2 color3;

varying vec2 vColor1;
varying vec2 vColor2;
varying vec2 vColor3;

void main() {
    vColor1 = color1;
    vColor2 = color2;
    vColor3 = color3;
    vPosition = (modelMatrix * vec4(position, 1.0)).xyz;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    #include <logdepthbuf_vertex>
}
`;

const fragmentShader = `
precision highp float;
precision highp int;

uniform float burnRadius;
varying vec3 vPosition;
varying vec2 vColor1;
varying vec2 vColor2;
varying vec2 vColor3;

#include <common>
#include <logdepthbuf_pars_fragment>

float logScale(float value) {
    float epsilon = 1e-6;  // Small constant to avoid log(0)
    float scaledValue = log(value + epsilon) / log(0.8 + epsilon);  // Log scale normalization
    return clamp(scaledValue, 0.0, 1.0);  // Keep the value in the range [0.0, 1.0]
}

void main() {
	#include <logdepthbuf_fragment>

    vec2 color1 = clamp(vColor1, 0.0, 1.0);
    vec2 color2 = clamp(vColor2, 0.0, 1.0);
    vec2 color3 = clamp(vColor3, 0.0, 1.0);

    gl_FragColor = vec4(color1, 0.0, 1.0);

    float burn = burnRadius;
    float cv1 = min(color1.r, color1.g) / burn;
    float cv2 = min(color2.r, color2.g) / burn;
    float cv3 = min(color3.r, color3.g) / burn;

    float cv = min(min(cv1, cv2), cv3);

    gl_FragColor = vec4(cv, cv, 0.0, 1.0);

    // float maxColor = logScale( (max(max(color.r, color.g), color.b) - 0.5) * 2.0 );
    // float maxColor = (max(max(color.r, color.g), color.b) - 0.5) * 2.0;
    // gl_FragColor = vec4(vec3(maxColor), 1.0);  // Apply grayscale with full opacity
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
// dep: moto.license
gapp.register("geo.toast", [], (root, exports) => {
    exports({
        create: createEdgeData
    });
});
