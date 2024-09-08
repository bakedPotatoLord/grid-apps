/** Copyright Stewart Allen <sa@grid.space> -- All Rights Reserved */

"use strict";

// toast/burnt-edge material and helpers

// modified EdgesGeometry to preserve edge/face relationships
function createEdgeData(geometry, thresholdAngle = 20) {

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
    const faceIdxs = new Array(numFaces); // [i1, i2] array per face

    function repl(arr, val) {
        if (arr.indexOf(val) >= 0) {
            return;
        }
        let io = arr.indexOf(0);
        if (io >= 0) {
            arr[io] = val;
        } else {
            console.log({ arr, io, val });
            throw "array full";
        }
    }

    for (let i = 0, faceId = 0; i < indexCount; i += 3, faceId++) {
        let idxs = faceIdxs[faceId] = [0,0,0,0];

        // points from face
        indexArr[0] = i;
        indexArr[1] = i + 1;
        indexArr[2] = i + 2;

        const { a, b, c } = _triangle;
        a.fromBufferAttribute(positionAttr, indexArr[0]);
        b.fromBufferAttribute(positionAttr, indexArr[1]);
        c.fromBufferAttribute(positionAttr, indexArr[2]);
        _triangle.getNormal(_normal);

        // create hashes for the edge from the vertices
        hashes[0] = `${Math.round(a.x * precision)},${Math.round(a.y * precision)},${Math.round(a.z * precision)}`;
        hashes[1] = `${Math.round(b.x * precision)},${Math.round(b.y * precision)},${Math.round(b.z * precision)}`;
        hashes[2] = `${Math.round(c.x * precision)},${Math.round(c.y * precision)},${Math.round(c.z * precision)}`;

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
            // if (reverseHash in edgeData && edgeData[reverseHash]) {
            if (adjacent) {
                // if we found a sibling edge add it into the vertex array if
                // it meets the angle threshold and delete the edge from the map.
                if (_normal.dot(edgeData[reverseHash].normal) <= thresholdDot) {
                    let lineIndex = vertices.length / 4 + 1;
                    vertices.push(v0.x, v0.y, v0.z, 0);
                    vertices.push(v1.x, v1.y, v1.z, 0);
                    repl(idxs, lineIndex);
                    // add line index to adjacent face
                    repl(faceIdxs[adjacent.faceId], lineIndex);
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

    // Create a BufferAttribute for the line indices
    const edgeIndexArr = faceIdxs.flat().toFloat32(); //new Float32Array(numFaces * 4);
    const edgeIndices = new THREE.DataTexture(edgeIndexArr, numFaces, 1, THREE.RGBAFormat, THREE.FloatType);
    const edgeLines = new THREE.DataTexture(vertices.toFloat32(), vertices.length/4, 1, THREE.RGBAFormat, THREE.FloatType);
    edgeIndices.needsUpdate = true;
    edgeIndices.minFilter = THREE.NearestFilter;
    edgeIndices.magFilter = THREE.NearestFilter;
    edgeLines.needsUpdate = true;
    edgeLines.minFilter = THREE.NearestFilter;
    edgeLines.magFilter = THREE.NearestFilter;

    const material = new ShaderMaterial({
        vertexShader: vertexShader,
        fragmentShader: fragmentShader,
        uniforms: {
            burnRadius: { value: 10.0 },
            edgeLines: { value: edgeLines },
            edgeIndices: { value: edgeIndices },
        }
    });

    const mgeo = geometry.clone();
    const mesh = new THREE.Mesh(mgeo, material);

    console.log({ numFaces, vertices, edgeLines, edgeIndices, faceIdxs: faceIdxs.flat() });

    return { vertices, material, mesh };
}

const vertexShader = `
uniform sampler2D edgeIndices;    // Texture containing edge indices (now a float texture)
uniform sampler2D edgeLines;      // Texture containing edge line endpoints (float texture)
varying vec3 vLineStart1, vLineEnd1, vLineStart2, vLineEnd2, vLineStart3, vLineEnd3;
varying float vLineCount;
varying vec3 vPosition;

vec3 getEdgeLine(int index) {
  // Ensure valid index (1-based to 0-based indexing)
  if (index < 1) return vec3(0.0);  // If no line, return dummy zero vector

  float texSize = float(textureSize(edgeLines, 0).x);  // Texture size (number of line segments)
  float u = float(index - 1) / texSize;  // Convert 1-based index to 0-based for texture lookup
  return texture2D(edgeLines, vec2(u, 0.0)).rgb;  // Sample the edgeLines texture
}

vec3 getEdgeIndices(int faceIndex) {
  float texSize = float(textureSize(edgeIndices, 0).x);  // Number of faces (texture width)
  float u = float(faceIndex) / texSize;  // Normalize face index for texture lookup
  return texture2D(edgeIndices, vec2(u, 0.0)).rgb;  // Sample edge indices as vec3
}

void main() {
  vLineCount = 0.0;

  // Calculate the face index using gl_VertexID
  int faceIndex = int(gl_VertexID) / 3;

  // Get edge indices for this face
  vec3 indices = getEdgeIndices(faceIndex);

  // Handle the first line
  if (indices.x > 0.0) {
    int index1 = int(indices.x);
    vLineStart1 = getEdgeLine(index1 * 2 - 1);  // Get the start of the first line
    vLineEnd1 = getEdgeLine(index1 * 2);        // Get the end of the first line
    vLineCount += 1.0;
  }

  // Handle the second line
  if (indices.y > 0.0) {
    int index2 = int(indices.y);
    vLineStart2 = getEdgeLine(index2 * 2 - 1);  // Get the start of the second line
    vLineEnd2 = getEdgeLine(index2 * 2);        // Get the end of the second line
    vLineCount += 1.0;
  }

  // Handle the third line
  if (indices.z > 0.0) {
    int index3 = int(indices.z);
    vLineStart3 = getEdgeLine(index3 * 2 - 1);  // Get the start of the third line
    vLineEnd3 = getEdgeLine(index3 * 2);        // Get the end of the third line
    vLineCount += 1.0;
  }

  // Store the transformed vertex position
  vPosition = (modelViewMatrix * vec4(position, 1.0)).xyz;

  // Pass position to the fragment shader
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

const fragmentShader = `
varying float vLineCount;  // Number of lines for this face
varying vec3 vLineStart1, vLineEnd1, vLineStart2, vLineEnd2, vLineStart3, vLineEnd3;
varying vec3 vPosition;  // Vertex position passed from the vertex shader
uniform float burnRadius;

float distanceToLine(vec3 point, vec3 start, vec3 end) {
  vec3 lineDir = normalize(end - start);
  vec3 v = point - start;
  float d = dot(v, lineDir);
  vec3 closestPoint = start + clamp(d, 0.0, length(end - start)) * lineDir;
  return length(closestPoint - point);
}

void main() {
  vec4 color = vec4(1.0, 1.0, 1.0, 1.0);  // Base color (white)

  float burnFactor = 0.0;

  // Calculate burn factor for the first line if present
  if (vLineCount >= 1.0) {
    float distToLine1 = distanceToLine(vPosition, vLineStart1, vLineEnd1);
    // burnFactor = max(burnFactor, smoothstep(burnRadius, 0.0, distToLine1));
    burnFactor += 0.25;
  }

  // Calculate burn factor for the second line if present
  if (vLineCount >= 2.0) {
    float distToLine2 = distanceToLine(vPosition, vLineStart2, vLineEnd2);
    // burnFactor = max(burnFactor, smoothstep(burnRadius, 0.0, distToLine2));
    burnFactor += 0.25;
  }

  // Calculate burn factor for the third line if present
  if (vLineCount >= 3.0) {
    float distToLine3 = distanceToLine(vPosition, vLineStart3, vLineEnd3);
    // burnFactor = max(burnFactor, smoothstep(burnRadius, 0.0, distToLine3));
    burnFactor += 0.25;
  }

  // Apply darkening effect based on burn factor
  vec4 burntColor = vec4(color.rgb * (1.0-burnFactor), color.a);  // Multiply only RGB, keep alpha the same

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
// dep: moto.license
gapp.register("geo.toast", [], (root, exports) => {
    exports({
        create: createEdgeData
    });
});
