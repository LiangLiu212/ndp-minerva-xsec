"""jobsub-lite skill library — generic jobsub_lite + RCDS helpers.

Vendored and generalized from genie-dev/jobsub-agent/lib (2026-07). Project
state (config, catalog, tarball cache, run records) lives in a per-project
`.jobsub/` directory resolved by lib.config.state_dir(); the skill directory
itself stays read-only so it can be copied between projects.
"""
