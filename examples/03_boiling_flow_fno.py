"""Neural-operator surrogate of nucleate boiling flow (BubbleML 2.0, single bubble, FC-72).

1. downloads a time window of five interface-resolved Flash-X simulations
   (heater temperatures 87-103 degC) directly from Hugging Face,
2. trains an autoregressive FNO on four of them,
3. rolls it out on the held-out heater temperature (94 degC) and compares the
   temperature/velocity/interface fields and the heater heat flux.

Requires:  pip install "boilflux[ml]"
Licence note: BubbleML 2.0 is CC BY-NC-SA 4.0. The data and any model trained
on it are for non-commercial research only and are not redistributed here.

Run:  python examples/03_boiling_flow_fno.py --out fno_run
"""

import argparse
import json
from pathlib import Path

from boilflux.operators.bubbleml import fetch_simulation
from boilflux.operators.train import (
    TrainConfig,
    evaluate_rollout,
    plot_rollout,
    save_report,
    summary,
    train_surrogate,
)

SUBSET = "SingleBubble-Saturated-FC72-2D"
TRAIN_T, TEST_T = (87, 91, 96, 103), 94

p = argparse.ArgumentParser()
p.add_argument("--data", default="bubbleml_cache")
p.add_argument("--out", default="fno_run")
p.add_argument("--epochs", type=int, default=25)
args = p.parse_args()

data, out = Path(args.data), Path(args.out)
files = {}
for T in (*TRAIN_T, TEST_T):
    f = data / f"sb_fc72_T{T}.npz"
    if not f.exists():
        print("downloading", SUBSET, T)
        fetch_simulation(SUBSET, T, f, frames=(600, 1200), stride=2, downsample=3)
    files[T] = f

model, info = train_surrogate([files[T] for T in TRAIN_T], TrainConfig(epochs=args.epochs))
out.mkdir(parents=True, exist_ok=True)
model.save(out / "fno_bubbleml_sb_fc72.pt")
ev = evaluate_rollout(model, files[TEST_T], steps=80)
plot_rollout(ev, out, steps_to_show=(0, 10, 25, 50, 79))
summaries = [summary(ev)]
save_report(info, summaries, out / "report.json")
print(json.dumps(summaries, indent=2))
