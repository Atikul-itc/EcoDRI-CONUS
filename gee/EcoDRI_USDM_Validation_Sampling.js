// ============================================================================
// EcoDRI USDM Validation Sample Generator
//
// Draws stratified random pixel samples at 80 evenly-spaced USDM Tuesday
// release dates across 2010-2024 for validation against the operational
// U.S. Drought Monitor classification.
//
// Output: one CSV per validation date, containing paired EcoDRI and USDM
// values, plus raw component values (VCI, TCI, SMCI, SPI) and Köppen zone
// for downstream ablation and stratified analyses.
// ============================================================================

var USDM = ee.ImageCollection('projects/sat-io/open-datasets/us-drought-monitor');
var ECODRI_COLL = 'projects/ee-atikul/assets/EcoDRI_v3_2_growing_season';

var STUDY_AREA = ee.Geometry.Rectangle([-125, 23, -65, 50]);
var SAMPLES_PER_CLASS = 200;   // per USDM category per date

// Köppen major group (see EcoDRI_Weekly_Composite.js)
var koppenShp = ee.FeatureCollection('projects/ee-atikul/assets/Aridity_Koppen');
var gridcodeImg = koppenShp.reduceToImage({
  properties: ['GRIDCODE'],
  reducer: ee.Reducer.first()
}).rename('gridcode');
var koppenMajor = gridcodeImg.divide(10).floor().toInt().rename('koppen');

// Generate 80 evenly-spaced Tuesday dates across 2010-2024 growing seasons
function generateValidationDates() {
  var dates = [];
  for (var y = 2010; y <= 2024; y++) {
    var d = new Date(Date.UTC(y, 3, 1));
    var dow = d.getUTCDay();
    d.setUTCDate(d.getUTCDate() + ((2 - dow + 7) % 7));   // first Tuesday
    while (d.getUTCMonth() <= 9) {
      dates.push(d.toISOString().slice(0, 10));
      d.setUTCDate(d.getUTCDate() + 8 * 7);   // 8-week step
    }
  }
  return dates;
}

function sampleOneDate(dateStr) {
  var d = ee.Date(dateStr);
  var win = 4;   // days

  var ecoImg = ee.Image(
    ee.ImageCollection(ECODRI_COLL)
      .filterDate(d.advance(-win,'day'), d.advance(win,'day'))
      .first()
  );

  var usdmImg = ee.Image(
    USDM.filterDate(d.advance(-win,'day'), d.advance(win,'day')).first()
  ).rename('usdm_cat').unmask(0);   // outside-polygon = None

  var stack = ecoImg
    .select(['drought_category','EcoDRI','VCI','TCI','SMCI','SPI'])
    .addBands(usdmImg)
    .addBands(koppenMajor);

  return stack.stratifiedSample({
    numPoints: SAMPLES_PER_CLASS,
    classBand: 'usdm_cat',
    region: STUDY_AREA,
    scale: 5000,
    seed: 42,
    dropNulls: true,
    tileScale: 8
  }).map(function(f) { return f.set('date', dateStr); });
}

var dates = generateValidationDates();
dates.forEach(function(dateStr) {
  var sample = sampleOneDate(dateStr);
  Export.table.toDrive({
    collection: sample,
    description: 'EcoDRI_USDM_paired_' + dateStr.replace(/-/g,'_'),
    folder: 'EcoDRI_USDM_Validation',
    fileFormat: 'CSV',
    selectors: ['date','koppen','usdm_cat','drought_category',
                'EcoDRI','VCI','TCI','SMCI','SPI']
  });
});

print('Submitted', dates.length, 'validation sampling tasks.');
