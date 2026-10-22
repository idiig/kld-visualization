// npm install --no-save playwright; npx playwright install chromium
// Serve _site, then: node tests/browser.mjs http://localhost:8080/
import fs from 'node:fs/promises';
import { spawn } from 'node:child_process';
const staticServer = process.env.KLD_SERVE ? spawn('python3', ['-m', 'http.server', '8080', '--bind', '127.0.0.1', '--directory', '_site'], { stdio: 'ignore' }) : null;
if (staticServer) await new Promise(resolve => setTimeout(resolve, 500));
const { chromium } = await import(process.env.PLAYWRIGHT_MODULE || 'playwright');
const browser = await chromium.launch({
  headless: true,
  executablePath: process.env.CHROMIUM_PATH || undefined,
  args: ['--no-sandbox', '--disable-dev-shm-usage'],
});
const page = await browser.newPage({ viewport: { width: 1440, height: 1100 }, acceptDownloads: true });
page.setDefaultTimeout(120000);
page.on('requestfailed', r => console.error('REQUEST FAILED:', r.url(), r.failure()));
page.on('pageerror', e => console.error('PAGE ERROR:', e.message));
page.on('console', msg => { if (msg.type() === 'error') console.error('CONSOLE:', msg.text()); });
try {
  await page.goto(process.argv[2] || 'http://localhost:8080/');
  let frame;
  for (let i = 0; i < 120; i++) {
    frame = page.frames().find(f => f.url().includes('index.html') && f !== page.mainFrame());
    if (!frame) frame = page.frames().find(f => f !== page.mainFrame());
    if (frame && await frame.locator('#render').count()) break;
    await page.waitForTimeout(1000);
  }
  if (!frame || !await frame.locator('#render').count()) throw new Error('App did not start');
  await frame.waitForFunction(() => window.Shiny?.shinyapp?.isConnected());
  console.log('Browser App ready:', frame.url());
  const choices = await frame.evaluate(() => Object.keys(document.querySelector('#input_file').selectize.options).length);
  if (choices !== 8) throw new Error(`Expected 8 datasets, found ${choices}`);
  await frame.evaluate(() => document.querySelector('#input_file').selectize.setValue('entry-vs-gender-kokin.tsv'));
  await frame.evaluate(() => document.querySelector('#badge_selection').selectize.setValue('frontier'));
  await frame.locator('#frontier_n').fill('3');
  await frame.locator('#render').click();
  await frame.waitForFunction(() => document.querySelector('#status')?.textContent.includes('描画完了'));
  const image = frame.locator('#preview img');
  await image.waitFor();
  await frame.waitForFunction(() => { const img = document.querySelector('#preview img'); return img?.complete && img.naturalWidth > 10; });
  const firstSource = await image.getAttribute('src');
  console.log('PASS: 8 datasets, parameter controls, PNG preview');
  const downloadEvent = page.waitForEvent('download');
  await frame.locator('#download_svg').click();
  const download = await downloadEvent;
  console.log('DOWNLOAD URL:', download.url());
  let svg;
  try { svg = await fs.readFile(await download.path(), 'utf8'); }
  catch (e) {
    console.error('DOWNLOAD FAILED:', await download.failure());
    const response = await frame.evaluate(async url => {
      const r = await fetch(url); return { status: r.status, text: await r.text() };
    }, download.url());
    console.error('DOWNLOAD RESPONSE:', response.status, response.text.slice(0, 1000));
    throw e;
  }
  if (!svg.includes('<svg') || svg.length < 1000) throw new Error('Invalid SVG download');
  console.log('PASS: SVG download', download.suggestedFilename(), svg.length);
  await frame.evaluate(() => document.querySelector('#kld_variant').selectize.setValue('obs'));
  await frame.locator('#render').click();
  await frame.waitForFunction(old => {
    const img = document.querySelector('#preview img');
    return img?.complete && img.naturalWidth > 10 && img.src !== old;
  }, firstSource);
  console.log('PASS: parameter change re-renders plot');
  await frame.locator('a[data-value="表"]').click();
  await frame.locator('#annotated_table tbody tr').first().waitFor();
  await frame.locator('a[data-value="プロット"]').click();
  await frame.evaluate(() => document.querySelector('#input_file').selectize.setValue('pair-entry-budatemajor-vs-gender-first4.tsv'));
  await frame.locator('#render').click();
  await frame.waitForFunction(() => document.querySelector('#status')?.textContent.includes('描画完了: pair-entry-budatemajor'));
  console.log('PASS: result table, pair dataset download and plot');
  await page.screenshot({ path: '/tmp/kld-browser-preview.png', fullPage: true });
} catch (e) {
  console.error('Frames:', page.frames().map(f => f.url()));
  for (const f of page.frames()) console.error((await f.locator('body').innerText().catch(() => '')).slice(-3000));
  await page.screenshot({ path: '/tmp/kld-browser-error.png', fullPage: true }).catch(() => {});
  throw e;
} finally {
  await browser.close();
  staticServer?.kill();
}
