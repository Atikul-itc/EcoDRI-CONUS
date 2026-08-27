// ============================================================================
// EcoDRI County Aggregation for Yield Validation
//
// Aggregates weekly EcoDRI to US county boundaries over two temporal windows:
//   - Full growing season (April - October)
//   - Crop-specific critical windows:
//       corn:    July 1 - August 15
//       soybean: July 20 - September 5
//       winter wheat: April 15 - June 5
// Output: one CSV per (year, crop, window).
// ============================================================================

var COUNTIES = ee.FeatureCollection('TIGER/2018/Counties');
var ECODRI_COLL = 'projects/ee-atikul/assets/EcoDRI_v3_2_growing_season';
var STUDY_AREA = ee.Geometry.Rectangle([-125, 23, -65, 50]);

var YEARS = ee.List.sequence(2010, 2024);
var CROP_WINDOWS = {
  corn:         {startMon: 7, startDay: 1,  endMon: 8, endDay: 15},
  soybean:      {startMon: 7, startDay: 20, endMon: 9, endDay: 5},
  winter_wheat: {startMon: 4, startDay: 15, endMon: 6, endDay: 5}
};

var conusCounties = COUNTIES.filterBounds(STUDY_AREA);

function aggregate(year, startMon, startDay, endMon, endDay, label) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, startMon, startDay);
  var end   = ee.Date.fromYMD(year, endMon, endDay).advance(1, 'day');

  var coll = ee.ImageCollection(ECODRI_COLL).filterDate(start, end);
  var meanImg = coll.select(['EcoDRI','VCI','TCI','SMCI','SPI']).mean()
    .rename(['EcoDRI_mean','VCI_mean','TCI_mean','SMCI_mean','SPI_mean']);
  var maxImg = coll.select('EcoDRI').max().rename('EcoDRI_max');
  var p90Img = coll.select('EcoDRI').reduce(ee.Reducer.percentile([90]))
    .rename('EcoDRI_p90');

  var stack = meanImg.addBands(maxImg).addBands(p90Img);

  return stack.reduceRegions({
    collection: conusCounties,
    reducer: ee.Reducer.mean(),
    scale: 5000
  }).map(function(f) {
    return f.set({year: year, window: label});
  });
}

// Full-season aggregation
YEARS.getInfo().forEach(function(y) {
  var summary = aggregate(y, 4, 1, 10, 31, 'full_season');
  Export.table.toDrive({
    collection: summary,
    description: 'EcoDRI_county_full_' + y,
    folder: 'EcoDRI_County_Aggregation',
    fileFormat: 'CSV',
    selectors: ['GEOID','STATEFP','NAME','year','window',
                'EcoDRI_mean','EcoDRI_max','EcoDRI_p90',
                'VCI_mean','TCI_mean','SMCI_mean','SPI_mean']
  });
});

// Crop-specific critical windows
YEARS.getInfo().forEach(function(y) {
  Object.keys(CROP_WINDOWS).forEach(function(crop) {
    var w = CROP_WINDOWS[crop];
    var summary = aggregate(y, w.startMon, w.startDay, w.endMon, w.endDay, crop);
    Export.table.toDrive({
      collection: summary,
      description: 'EcoDRI_county_' + crop + '_' + y,
      folder: 'EcoDRI_County_Aggregation',
      fileFormat: 'CSV',
      selectors: ['GEOID','STATEFP','NAME','year','window',
                  'EcoDRI_mean','EcoDRI_max','EcoDRI_p90']
    });
  });
});
