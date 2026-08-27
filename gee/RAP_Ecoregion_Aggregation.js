// ============================================================================
// RAP Herbaceous Biomass Ecoregion Aggregation
//
// Reads the Rangeland Analysis Platform 16-day NPP product and aggregates
// peak-season (May-September) herbaceous biomass (afgNPP + pfgNPP) to EPA
// Level III ecoregions for the ecological validation.
// ============================================================================

var STUDY_AREA = ee.Geometry.Rectangle([-125, 23, -65, 50]);
var RANGELAND_FRACTION_THRESHOLD = 0.30;
var LBS_PER_ACRE_TO_KG_PER_HA = 1.121;

var rapNPP = ee.ImageCollection(
  'projects/rangeland-analysis-platform/npp-partitioned-16day-v3');

var ecoregionsL3 = ee.FeatureCollection('EPA/Ecoregions/2013/L3').filterBounds(STUDY_AREA);
var nlcd = ee.Image('USGS/NLCD_RELEASES/2019_REL/NLCD/2019').select('landcover');
var rangelandMask = nlcd.eq(52).or(nlcd.eq(71)).rename('rangeland');

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

var YEARS = ee.List.sequence(2010, 2024);

function aggregateYear(year) {
  year = ee.Number(year);
  var start = ee.Date.fromYMD(year, 5, 1);
  var end   = ee.Date.fromYMD(year, 9, 30).advance(1, 'day');
  var coll = rapNPP.filterDate(start, end);

  var afg = coll.select('afgNPP').sum().multiply(LBS_PER_ACRE_TO_KG_PER_HA);
  var pfg = coll.select('pfgNPP').sum().multiply(LBS_PER_ACRE_TO_KG_PER_HA);
  var shr = coll.select('shrNPP').sum().multiply(LBS_PER_ACRE_TO_KG_PER_HA);
  var herb = afg.add(pfg).rename('NPP_herbaceous');
  var range = herb.add(shr).rename('NPP_rangeland');

  var stack = herb.addBands(range)
    .addBands(afg.rename('NPP_annual_herb'))
    .addBands(pfg.rename('NPP_perennial_herb'))
    .updateMask(rangelandMask);

  return stack.reduceRegions({
    collection: rangelandEcoregions,
    reducer: ee.Reducer.mean(),
    scale: 1000
  }).map(function(f) { return f.set('year', year); });
}

YEARS.getInfo().forEach(function(y) {
  var summary = aggregateYear(y);
  Export.table.toDrive({
    collection: summary,
    description: 'RAP_ecoregion_peak_' + y,
    folder: 'RAP_Ecoregion_Aggregation',
    fileFormat: 'CSV',
    selectors: ['na_l3code','na_l3name','na_l1name','year','rangeland_frac',
                'NPP_herbaceous','NPP_rangeland',
                'NPP_annual_herb','NPP_perennial_herb']
  });
});
