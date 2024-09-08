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

    const material = new ShaderMaterial({
        transparent: true,
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
        uniforms: {
            burnRadius: { value: 5.0 },
        }
    });

    const colorAttribute = positionAttr.clone();
    const colors = colorAttribute.array;
    const zero = -0.0;
    const R = 1.0, G = 0.0, B = 0.0;
    for (let i = 0; i < colors.length; ) {
        colors[i++] = R * 1.0;
        colors[i++] = G * 1.0;
        colors[i++] = B * zero;

        colors[i++] = R * zero;
        colors[i++] = G * 1.0;
        colors[i++] = B * 1.0;

        colors[i++] = R * 1.0;
        colors[i++] = G * zero;
        colors[i++] = B * 1.0;
    }

    const mgeo = geometry.clone();
    mgeo.setAttribute('color', colorAttribute);
    const mesh = new THREE.Mesh(mgeo, material);

    console.log({ numFaces, faceRecords, pointToEdge, edges: vertices.group(4) });

    return { vertices, material, mesh };
}

const vertexShader = `
attribute vec3 color;
varying vec3 vPosition;
varying vec3 vColor;

void main() {
    vColor = color;
    vPosition = (modelMatrix * vec4(position, 1.0)).xyz;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

const fragmentShader = `
uniform float burnRadius;         // Proximity radius for burn effect
varying vec3 vPosition;           // Vertex position passed from the vertex shader
varying vec3 vColor;              // Interpolated color from vertex shader

float logScale(float value) {
    float epsilon = 1e-6;  // Small constant to avoid log(0)
    float scaledValue = log(value + epsilon) / log(0.8 + epsilon);  // Log scale normalization
    return clamp(scaledValue, 0.0, 1.0);  // Keep the value in the range [0.0, 1.0]
}

void main() {
    vec3 color = clamp(vColor, 0.0, 1.0);
    gl_FragColor = vec4(color, 1.0);  // stock color gradient

    // float maxColor = logScale( (max(max(color.r, color.g), color.b) - 0.5) * 2.0 );
    float maxColor = (max(max(color.r, color.g), color.b) - 0.5) * 2.0;
    gl_FragColor = vec4(vec3(maxColor), 1.0);  // Apply grayscale with full opacity
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
