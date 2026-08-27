// ============================================================================
// EcoDRI Ecoregion Aggregation for Ecological Validation
//
// Aggregates weekly EcoDRI to EPA Level III ecoregions over the peak
// growing season (May - September), restricted to rangeland-dominated
// ecoregions (NLCD rangeland fraction >= 30%) and masked to natural
// rangeland pixels only.
//
// Produces the CDI side of the RAP/RCMAP ecological validation.
// ============================================================================

var ECODRI_COLL = 'projects/ee-atikul/assets/EcoDRI_v3_2_growing_season';
var STUDY_AREA = ee.Geometry.Rectangle([-125, 23, -65, 50]);
var RANGELAND_FRACTION_THRESHOLD = 0.30;

var ecoregionsL3 = ee.FeatureCollection('EPA/Ecoregions/2013/L3').filterBounds(STUDY_AREA);
var nlcd = ee.Image('USGS/NLCD_RELEASES/2019_REL/NLCD/2019').select('landcover');
var rangelandMask = nlcd.eq(52).or(nlcd.eq(71)).rename('rangeland');

// Filter to rangeland-dominated ecoregions
var ecoregionsWithFraction = ecoregionsL3.map(function(eco) {
  var frac = rangelandMask.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: eco.geometry(),
    scale: 1000,
    maxPixels: 1e10,
    bestEffort: true
  }).get('rangeland');
  return eco.set('rangeland_frac', frac);
});

var rangelandEcoregions = ecoregionsWithFraction
  .filter(ee.Filter.gte('rangeland_frac', RANGELAND_FRACTION_THRESHOLD));

print('Rangeland-dominated ecoregions:', rangelandEcoregions.size());

var YEARS = ee.List.sequence(2010, 2024);

function aggregatePeakSeason(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 5, 1);
  var end   = ee.Date.fromYMD(year, 9, 30).advance(1, 'day');

  var coll = ee.ImageCollection(ECODRI_COLL).filterDate(start, end);

  var meanImg = coll.select(['EcoDRI','VCI','TCI','SMCI','SPI']).mean()
    .updateMask(rangelandMask)
    .rename(['EcoDRI_mean','VCI_mean','TCI_mean','SMCI_mean','SPI_mean']);
  var maxImg = coll.select('EcoDRI').max()
    .updateMask(rangelandMask).rename('EcoDRI_max');

  return meanImg.addBands(maxImg).reduceRegions({
    collection: rangelandEcoregions,
    reducer: ee.Reducer.mean(),
    scale: 5000
  }).map(function(f) { return f.set('year', year); });
}

YEARS.getInfo().forEach(function(y) {
  var summary = aggregatePeakSeason(y);
  Export.table.toDrive({
    collection: summary,
    description: 'EcoDRI_ecoregion_peak_' + y,
    folder: 'EcoDRI_Ecoregion_Aggregation',
    fileFormat: 'CSV',
    selectors: ['na_l3code','na_l3name','na_l2code','na_l1code','na_l1name',
                'year','rangeland_frac',
                'EcoDRI_mean','EcoDRI_max',
                'VCI_mean','TCI_mean','SMCI_mean','SPI_mean']
  });
});
