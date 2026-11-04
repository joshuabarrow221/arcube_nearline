# Offline GIF encoder

`gifenc-1.0.3.js` is gifenc 1.0.3 from the npm package maintained at
https://github.com/mattdesl/gifenc, under the MIT license.

The upstream CommonJS `dist/gifenc.js` is wrapped in a local scope exposing
`globalThis.GIFEnc`; its external source-map directive is removed. The license
is embedded in the bundle so standalone HTML exports retain the notice.
No network access or package manager is needed at viewer runtime.

Upstream npm tarball SHA-1: `3a1a97fa5ab70ebefd000350b9aac8313b30b3aa`.
