// ============================================================================
// EcoDRI Batch Export
//
// Generates weekly EcoDRI composites (April through October) for all years
// in the study period and submits export tasks to Google Drive OR to a
// persistent Earth Engine ImageCollection asset.
//
// Runtime: submits ~720 export tasks (24 years × 30 weeks/year). Google
// Earth Engine will queue them and process serially over ~24 hours.
// ============================================================================

var composite = require('users/atikul/EcoDRI:EcoDRI_Weekly_Composite');

var START_YEAR = 2001;
var END_YEAR   = 2024;
var STUDY_AREA = ee.Geometry.Rectangle([-125, 23, -65, 50]);
var OUTPUT_ASSET_COLLECTION = 'projects/ee-atikul/assets/EcoDRI_v3_2_growing_season';

// Generate one weekly composite per Wednesday from April through October
function generateDates(year) {
  var start = ee.Date.fromYMD(year, 4, 1);
  var dates = ee.List.sequence(0, 30).map(function(i) {
    return start.advance(ee.Number(i).multiply(7), 'day').millis();
  });
  return dates.getInfo().map(function(ms) {
    return new Date(ms).toISOString().slice(0, 10);
  });
}

for (var y = START_YEAR; y <= END_YEAR; y++) {
  var dates = generateDates(y);
  dates.forEach(function(dateStr) {
    var img = composite.generateWeekly(dateStr);
    var assetId = OUTPUT_ASSET_COLLECTION + '/EcoDRI_' + dateStr.replace(/-/g,'_');
    Export.image.toAsset({
      image: img,
      description: 'EcoDRI_' + dateStr,
      assetId: assetId,
      region: STUDY_AREA,
      scale: 5000,
      crs: 'EPSG:4326',
      maxPixels: 1e13
    });
  });
}

print('Submitted', (END_YEAR - START_YEAR + 1) * 31, 'export tasks.');
