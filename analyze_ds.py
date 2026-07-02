import json

with open('libs/rcad-algorithms/pipeline_dumps/rcad_bfuse_simple_A1_s03_after_BuildResultWire_ds.json', 'r') as f:
    ds = json.load(f)

print('=== ICs ===')
for ic in ds['ds']['intersection_curves']:
    ci = ic['ci']
    sv = ic['sv']
    ev = ic['ev']
    print(f'IC {ci}: sv={sv} ev={ev} nPBs={ic["n_pave_blocks"]}')
    print(f'  curve: {str(ic["curve"])[:80]}')

print()
print('=== Edge vertex pairs ===')
for e in ds['ds']['edges']:
    sv = e['sv']
    ev = e['ev']
    print(f'  ei={e["ei"]} v{sv} -> v{ev}  {str(e["curve"])[:60]}')

print()
print('=== Face curves_sc ===')
for f in ds['ds']['faces']:
    print(f'  fi={f["fi"]} nPBsSc={f["nPBsSc"]} nCurvesSc={f["nCurvesSc"]} boundary={f["boundary_edges"]}')
