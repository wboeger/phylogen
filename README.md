# phylogen

Build a phylogenetic tree from any species list — free, open, general-purpose.

Fetches sequences from GenBank (NCBI), aligns them (MAFFT via usegalaxy.eu),
trims the alignment (trimAl via usegalaxy.eu), and produces an instant local
Neighbor-Joining preview tree, with an optional bootstrapped Maximum-Likelihood
tree (RAxML-NG on usegalaxy.eu) on request. Works for any taxon NCBI has
sequence data for — not limited to one clade.

Adapted from the phylogeny module of
[AI_morpho / GyroMorpho](https://github.com/wboeger/AI_morpho) (a Gyrodactylidae
morphometrics tool), generalized and stripped of everything specific to that
app's Specimen/Project morphometrics database.

## How it works

- No login. A signed browser cookie remembers which jobs are yours (for the
  homepage's "recent jobs" list); every job also gets a shareable URL.
- No shared server secrets. Every job supplies its own NCBI contact email
  (+ optional API key) and its own usegalaxy.eu API key. Credentials are held
  in an in-process memory cache for the lifetime of the pipeline call and are
  **never written to disk or the database**.
- Ephemeral storage only — no persistent volume. Job history and files live
  for the life of the running container; a redeploy wipes them. The real
  deliverable is the ZIP a finished job lets you download to your own machine.

## Local development

```bash
python3 -m venv venv
./venv/bin/pip install -r requirements.txt
./venv/bin/python run.py   # http://127.0.0.1:5001
```

## Deploy (Railway)

Nixpacks auto-detects Python; `Procfile` / `railway.json` define the start
command (`gunicorn`, 1 worker / 8 threads). No volume needed. Environment
variables (all optional — every credential the pipeline needs comes from the
user, per job):

| Variable | Purpose |
|---|---|
| `SECRET_KEY` | Flask session-cookie signing (set to a long random string) |
| `DATA_DIR` | Override the ephemeral working directory (defaults to `./data`) |
| `GALAXY_BASE_URL` | Galaxy server (defaults to `https://usegalaxy.eu`) |
| `GALAXY_MAFFT_TOOL_ID` / `GALAXY_TRIMAL_TOOL_ID` / `GALAXY_RAXMLNG_TOOL_ID` | Pin Galaxy tool versions if the default server updates them |

Generate a public domain: Service → Settings → Networking → Generate Domain.
