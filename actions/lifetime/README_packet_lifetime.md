# Packet-based electron lifetime

The implementation follows Ahmad's packet-analysis workflow. The monitoring
profile extracts detector-wide charge objects from native packets
and pools them in fixed six-hour windows before fitting. Accepted results and
finite failed-check candidates appear on the main plot. See the
[monitoring README](README.md) for interpretation and
[methods reference](../../docs/lifetime-methods.md#packets) for gates and errors.

## Inputs and calibration

The extractor needs a native packet HDF5 and its matching FLOW file. It subtracts
static per-channel pedestals, groups successive packets from the same physical
channel, sums signed charge, and uses the last packet's EXT-relative time for
each object. The wrapper uses the 28-tick grouping rule and EXT I/O group configured
for the installation. These timing/grouping choices remain provisional.

The configured monitor uses the exact static pedestal JSON named in FLOW and
checks it against linked charge/ADC pairs. Missing calibration channels remain
missing and affect coverage. The standalone command also supports recovering
pedestals directly from FLOW (`--ped-source flow`) or explicitly supplied legacy
panel pedestals (`--ped-source panel`). These are distinct calibration policies;
do not combine their objects silently.

## Standalone extraction

Install [requirements.txt](requirements.txt). To extract reusable objects using
the same static-calibration policy as the monitor:

```bash
python packet_lifetime.py /data/run/packet.h5 \
  --flow-file /data/run/packet.FLOW.hdf5 \
  --ped-source calibration --ped /calibration/pedestal.json \
  --detector-wide --charge-mode adc-minus-ped \
  --successive-min-ticks 28 --successive-max-ticks 28 \
  --ext-io 6 --cohort verified-stable-run \
  --sample-timestamp "$ACQUISITION_START_ISO" \
  --sample-directory /output/samples --outdir /output/extraction --extract-only
```

Use the actual acquisition timestamp, matching calibration and verified EXT
settings. Then run a declared-window study with
`packet_pooling.py --samples /output/samples --hours 6 --bootstrap 100 --output /output/pools.json`.
Use each command's `--help` for alternatives. The standalone shard store identifies
files by basename; keep unrelated runs with colliding names separate. The
configured monitor instead uses full relative acquisition identities.

## Fit status and resources

Charge-slice fits, bootstrap stability, drift coverage, attenuation quality and
per-file pedestal coverage all contribute to status. Pooling increases statistics
but does not repair missing calibration or a wrong charge model. Fixed gates and
window membership are retained regardless of agreement with other methods.
Failed candidates keep their status when displayed; missing errors are not zero.

Complete native timing arrays are loaded into memory, so worker RAM must suit the
largest file; `--chunk` is not a total-memory bound. Extraction diagnostics and
compact hashed objects remain available for audit/retry. A packet failure does
not prevent available PRM, gas and track results from updating.

The track comparison curves overlap: cosmic-enriched tracks are the
`n_ext_trigs == 0` subset of All tracks. They share selected segments and are not
two independent validations of a packet estimate. See the
[track definitions](README.md#track-selections-and-trigger-categories).
