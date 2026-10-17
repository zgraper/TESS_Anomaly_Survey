// Optional Node acceptance test of the actual embedded browser function.
// Usage: node tests/progress_race.cjs [path/to/installed/web.py]
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const source = fs.readFileSync(process.argv[2] || path.join(__dirname, '../src/tess_anomaly_survey/web.py'), 'utf8');
const start = source.indexOf('    async function pollProgressOnce()');
const end = source.indexOf('    function updateProgressUI(', start);
assert(start >= 0 && end > start);
let release;
let updates = 0;
let reloads = 0;
const context = vm.createContext({
  AppState: {activeRun: {run_id: 'old', status: 'running'}},
  fetch: () => new Promise(resolve => { release = resolve; }),
  updateProgressUI: () => { updates++; },
  loadStats: () => { reloads++; },
  fetchStars: () => { reloads++; },
  clearInterval: () => {},
  console,
});
vm.runInContext(source.slice(start, end), context);
(async () => {
  const pending = context.pollProgressOnce();
  context.AppState.activeRun = {run_id: 'new', status: 'completed'};
  release({json: async () => ({run: {run_id: 'old', status: 'completed'}})});
  await pending;
  assert.equal(context.AppState.activeRun.run_id, 'new');
  assert.equal(updates, 0);
  assert.equal(reloads, 0);
  const current = context.pollProgressOnce();
  release({json: async () => ({run: {run_id: 'new', status: 'completed'}})});
  await current;
  assert.equal(updates, 1);
  assert.equal(reloads, 2);
  console.log('PASS: stale progress cannot replace selected run or reload its catalog; current response updates normally');
})().catch(error => { console.error(error); process.exitCode = 1; });
