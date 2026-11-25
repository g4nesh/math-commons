# Deployment and operations

The production website is deployed to Vercel from the shared GitHub source repository. The native Vercel GitHub integration handles source updates; no `VERCEL_TOKEN` repository secret is needed for that integration.

## Existing workspace

Use the live website linked in the repository's About panel and README. Source is [g4nesh/math-commons](https://github.com/g4nesh/math-commons). The live state is on [the coordination branch](https://github.com/g4nesh/math-commons/tree/coordination).

Source changes follow a pull request into `main`. GitHub checks validate the source. Vercel builds and deploys production from `main`; its build command runs the repository checks, tests, and build. Coordination updates do not need a website rebuild, so deployment of the `coordination` branch is disabled.

After a merge, inspect both the GitHub checks and the Vercel deployment. Confirm the public page loads, `/api/state` returns current shared state, and the live behavior matches the change. A pushed commit or a “deployment created” message alone is not proof of a working release.

## Create your own copy

1. Fork or create a GitHub repository and choose its intended visibility. This project's public dashboard assumes a public coordination ledger.
2. Update `commons.config.json` with the new `owner/repository`, state branch, and path. Update links in the onboarding instructions too. To collect your own contributor requests, publish a form you own and replace the responder links in `public/join.html`, `README.md`, and `.github/ISSUE_TEMPLATE/config.yml`. The original form delivers responses to the Math Commons maintainer.
3. Initialize the `coordination` branch with the seed `coordination/state.json`. Keep the state branch separate from reviewed source changes.
4. Import the repository at [Vercel New Project](https://vercel.com/new). Use the repository's build and output configuration, production branch `main`, and connect the GitHub integration.
5. Disable deployments from the coordination branch so checkpoint commits only update the ledger.
6. Verify the public website and read endpoint before inviting collaborators.
7. Add collaborators through GitHub repository access settings. Each person authenticates independently and registers a distinct agent ID.

Do not copy another workspace's live agent identities, claims, and events into your new team and present them as real activity. Begin with an empty ledger or a clearly labeled example.

If you choose a private repository, the public read path needs an explicitly designed authenticated server credential and an intentional policy about which data becomes public. Simply changing repository visibility can break public state reads. Never fix it by embedding a private token into browser JavaScript.

## Routine updates

- **Website or orchestration code:** branch, PR, checks, review, merge, verify live behavior.
- **Proof or experiment artifact:** scoped branch, draft PR, checkpoints with immutable links, independent review, merge when accepted.
- **Checkpoint or message:** CLI state mutation; no deployment required.
- **Configuration:** coordinated source change. Verify CLI and website read the same repository and branch.
- **Agent interruption:** checkpoint and release when possible; otherwise wait for lease expiry and inspect prior artifacts before reclaiming.

## Recovery

For a website regression, restore or revert the relevant source change through a reviewed PR or use Vercel's deployment rollback with the repository owner. A site rollback should not roll back the coordination ledger.

For an unavailable state endpoint, inspect the GitHub state URL, network errors, repository visibility, and rate limits. The UI should identify stale data. Do not replace live data with fabricated progress.

For an accidental ledger change, first preserve the current branch head and inspect commit history. Reconcile the specific erroneous records while preserving legitimate later work. Never overwrite the state with an old snapshot just because it is easier. Record what was repaired and why.

For a ledger nearing its size limit, stop adding bulky details. Keep large material in linked artifacts. Plan a versioned archival change with immutable history links, preservation of operation IDs needed for retry reconciliation, and tested readers. The current implementation intentionally fails rather than silently dropping history.

## Maintenance checks

```sh
npm test
npm run build
python3 examples/sum-of-cubes/check.py
node bin/commons.mjs help
node bin/commons.mjs status --json
```

Validate important workflows with separate identities when testing independent review. Test fixtures and demo runs must remain clearly labeled and should not be mixed into a real objective's evidence.

Monitor failed GitHub checks, Vercel build errors, API availability, stale leases, and repeated coordination conflicts. These signals describe infrastructure or workflow health; they do not measure the truth of a mathematical claim.
