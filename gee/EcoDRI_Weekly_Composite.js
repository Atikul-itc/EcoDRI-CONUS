// ============================================================================
// EcoDRI Weekly Composite Generator
//
// Produces weekly Ecological Drought Index (EcoDRI) composites over CONUS
// following the methodology described in Islam et al. (2024). Each weekly
// composite integrates:
//   - VCI  from MOD09A1 NDVI (500 m, 8-day)
//   - TCI  from MOD11A2 LST  (1 km, 8-day)
//   - SMCI from GLDAS + SMAP soil moisture
//   - SPI  from GRIDMET 90-day Standardized Precipitation Index
// The four components are combined into a single EcoDRI scalar using
// Köppen-Geiger climate-zone-specific weights, then discretized into six
// categorical drought classes aligned with the U.S. Drought Monitor.
//
// This script is invoked by EcoDRI_Batch_Export.js. Alone it renders one
// week for inspection.
//
// Reproducibility: all input datasets and version tags are pinned in the
// section headers below. Modify only the START_DATE and END_DATE variables
// for date-range changes.
// ============================================================================

// ---- CONFIGURATION ---------------------------------------------------------
var STUDY_AREA = ee.Geometry.Rectangle([-125, 23, -65, 50]);
var OUTPUT_SCALE = 5000;   // 5 km, Plate Carrée
var OUTPUT_CRS   = 'EPSG:4326';

// Week-of-year percentile climatology assets (pre-computed 2001-2024)
var CLIM_NDVI = 'projects/ee-atikul/assets/climatology_p5p95_NDVI_woy';
var CLIM_LST  = 'projects/ee-atikul/assets/climatology_p5p95_LST_woy';
var CLIM_SM   = 'projects/ee-atikul/assets/climatology_p5p95_SM_woy';

// Köppen-Geiger classification (user-provided asset)
var KOPPEN_SHP = 'projects/ee-atikul/assets/Aridity_Koppen';

// Recalibrated categorical thresholds (empirically calibrated against USDM;
// see Islam et al. 2024, Table 2)
var THRESHOLDS = [1.50, 1.85, 2.40, 2.90, 3.15];

// Epsilon guards on component ranges to prevent division-by-zero artifacts
var EPS_NDVI = 0.02;
var EPS_LST  = 0.5;
var EPS_SM   = 0.02;

// ---- INPUT COLLECTIONS -----------------------------------------------------
var NDVI_COLL = ee.ImageCollection('MODIS/061/MOD09A1')
  .select(['sur_refl_b01','sur_refl_b02','StateQA']);

var LST_COLL = ee.ImageCollection('MODIS/061/MOD11A2')
  .select(['LST_Day_1km','QC_Day']);

var GLDAS = ee.ImageCollection('NASA/GLDAS/V021/NOAH/G025/T3H')
  .select('SoilMoi0_10cm_inst');

var SMAP = ee.ImageCollection('NASA_USDA/HSL/SMAP10KM_soil_moisture')
  .select('ssm');

var GRIDMET_DROUGHT = ee.ImageCollection('GRIDMET/DROUGHT')
  .select('spi90d');

// ---- KÖPPEN CLIMATE-ZONE WEIGHTS -------------------------------------------
// Two-digit GRIDCODE: first digit is the major group.
// A=Tropical, B=Arid, C=Temperate, D=Continental, E=Polar
var koppenShp = ee.FeatureCollection(KOPPEN_SHP);
var gridcodeImg = koppenShp.reduceToImage({
  properties: ['GRIDCODE'],
  reducer: ee.Reducer.first()
}).rename('gridcode');
var koppenMajor = gridcodeImg.divide(10).floor().toInt().rename('koppen');

// Component weights by Köppen major group [VCI, TCI, SMCI, SPI]
// (Table 1 of Islam et al. 2024)
function koppenWeights(k) {
  var wVCI = ee.Image(0.25);
  var wTCI = ee.Image(0.25);
  var wSMCI = ee.Image(0.25);
  var wSPI = ee.Image(0.25);
  // A tropical
  wVCI  = wVCI.where(k.eq(1), 0.20);
  wTCI  = wTCI.where(k.eq(1), 0.20);
  wSMCI = wSMCI.where(k.eq(1), 0.20);
  wSPI  = wSPI.where(k.eq(1), 0.40);
  // B arid
  wVCI  = wVCI.where(k.eq(2), 0.15);
  wTCI  = wTCI.where(k.eq(2), 0.40);
  wSMCI = wSMCI.where(k.eq(2), 0.15);
  wSPI  = wSPI.where(k.eq(2), 0.30);
  // C temperate
  wVCI  = wVCI.where(k.eq(3), 0.20);
  wTCI  = wTCI.where(k.eq(3), 0.25);
  wSMCI = wSMCI.where(k.eq(3), 0.20);
  wSPI  = wSPI.where(k.eq(3), 0.35);
  // D continental
  wVCI  = wVCI.where(k.eq(4), 0.25);
  wTCI  = wTCI.where(k.eq(4), 0.25);
  wSMCI = wSMCI.where(k.eq(4), 0.20);
  wSPI  = wSPI.where(k.eq(4), 0.30);
  // E polar
  wVCI  = wVCI.where(k.eq(5), 0.10);
  wTCI  = wTCI.where(k.eq(5), 0.30);
  wSMCI = wSMCI.where(k.eq(5), 0.30);
  wSPI  = wSPI.where(k.eq(5), 0.30);
  return {vci: wVCI, tci: wTCI, smci: wSMCI, spi: wSPI};
}

// ---- COMPONENT COMPUTATION -------------------------------------------------
function maskMOD09A1(img) {
  var qa = img.select('StateQA');
  var clear = qa.bitwiseAnd(3).eq(0)
    .and(qa.bitwiseAnd(1 << 2).eq(0))
    .and(qa.bitwiseAnd(3 << 6).eq(0));
  return img.updateMask(clear);
}

function maskMOD11A2(img) {
  var qc = img.select('QC_Day');
  var good = qc.bitwiseAnd(3).lte(1);
  return img.updateMask(good);
}

function computeNDVI(img) {
  return img.normalizedDifference(['sur_refl_b02','sur_refl_b01']).rename('NDVI');
}

// VCI, TCI, SMCI from percentile climatology
function computeVCI(ndvi, woy) {
  var clim = ee.Image(CLIM_NDVI).select(['p5_' + woy, 'p95_' + woy])
    .rename(['p5','p95']);
  var range = clim.select('p95').subtract(clim.select('p5'))
    .max(ee.Image(EPS_NDVI));
  return ndvi.subtract(clim.select('p5')).divide(range)
    .multiply(100).clamp(0, 100).rename('VCI');
}

function computeTCI(lst, woy) {
  var clim = ee.Image(CLIM_LST).select(['p5_' + woy, 'p95_' + woy])
    .rename(['p5','p95']);
  var range = clim.select('p95').subtract(clim.select('p5'))
    .max(ee.Image(EPS_LST));
  return clim.select('p95').subtract(lst).divide(range)
    .multiply(100).clamp(0, 100).rename('TCI');
}

function computeSMCI(sm, woy) {
  var clim = ee.Image(CLIM_SM).select(['p5_' + woy, 'p95_' + woy])
    .rename(['p5','p95']);
  var range = clim.select('p95').subtract(clim.select('p5'))
    .max(ee.Image(EPS_SM));
  return sm.subtract(clim.select('p5')).divide(range)
    .multiply(100).clamp(0, 100).rename('SMCI');
}

// ---- WEEKLY COMPOSITE ------------------------------------------------------
function generateWeekly(targetDate) {
  var d = ee.Date(targetDate);
  var woy = d.getRelative('week', 'year').add(1).format('%02d').getInfo();

  // NDVI: nearest 8-day MODIS composite in ±4 day window
  var ndviWin = NDVI_COLL.filterDate(d.advance(-4,'day'), d.advance(4,'day'))
    .map(maskMOD09A1).map(computeNDVI);
  var ndvi = ndviWin.mean();

  // LST: same strategy
  var lstWin = LST_COLL.filterDate(d.advance(-4,'day'), d.advance(4,'day'))
    .map(maskMOD11A2);
  var lst = lstWin.select('LST_Day_1km').mean()
    .multiply(0.02);   // scale factor

  // Soil moisture: SMAP preferred, GLDAS fallback (±7 day window)
  var sm7 = d.advance(-7,'day');
  var e7  = d.advance(7,'day');
  var smapImg = SMAP.filterDate(sm7, e7).mean();
  var gldasImg = GLDAS.filterDate(sm7, e7).mean()
    .divide(50);   // ~saturation normalization
  var sm = ee.Image(ee.Algorithms.If(
    d.millis().gte(ee.Date('2015-04-01').millis()),
    smapImg.unmask(gldasImg),
    gldasImg
  ));

  // SPI: nearest daily GRIDMET DROUGHT observation
  var spi = GRIDMET_DROUGHT.filterDate(sm7, e7).mean().rename('SPI');

  // Convert to condition indices
  var vci  = computeVCI(ndvi, woy);
  var tci  = computeTCI(lst, woy);
  var smci = computeSMCI(sm, woy);

  // Drought-direction normalization D_i ∈ [0, 1]
  var dVCI  = ee.Image(1).subtract(vci.divide(100));
  var dTCI  = ee.Image(1).subtract(tci.divide(100));
  var dSMCI = ee.Image(1).subtract(smci.divide(100));
  var dSPI  = spi.multiply(-1).add(3).divide(6).clamp(0, 1);

  // Köppen-weighted aggregation with per-pixel weight renormalization
  // (handles pixels where one or more components are masked)
  var w = koppenWeights(koppenMajor);
  var mVCI  = dVCI.mask();
  var mTCI  = dTCI.mask();
  var mSMCI = dSMCI.mask();
  var mSPI  = dSPI.mask();
  var sumW = w.vci.multiply(mVCI)
    .add(w.tci.multiply(mTCI))
    .add(w.smci.multiply(mSMCI))
    .add(w.spi.multiply(mSPI));
  var num = dVCI.unmask(0).multiply(w.vci).multiply(mVCI)
    .add(dTCI.unmask(0).multiply(w.tci).multiply(mTCI))
    .add(dSMCI.unmask(0).multiply(w.smci).multiply(mSMCI))
    .add(dSPI.unmask(0).multiply(w.spi).multiply(mSPI));
  var ecoDRI = num.divide(sumW).multiply(4).rename('EcoDRI');

  // 3x3 focal-mean gap-fill
  ecoDRI = ecoDRI.focalMean(2.5, 'circle', 'pixels', 1).blend(ecoDRI);

  // Categorical classification
  var cat = ee.Image(0).rename('drought_category');
  for (var i = 0; i < THRESHOLDS.length; i++) {
    cat = cat.where(ecoDRI.gte(THRESHOLDS[i]), i + 1);
  }

  // Auxiliary VHI band (Kogan 2001; not used in EcoDRI, provided for reference)
  var vhi = vci.multiply(0.5).add(tci.multiply(0.5)).rename('VHI');

  return ecoDRI
    .addBands(cat)
    .addBands(vci)
    .addBands(tci)
    .addBands(smci)
    .addBands(vhi)
    .addBands(spi)
    .set('system:time_start', d.millis())
    .set('week_of_year', woy);
}

exports.generateWeekly = generateWeekly;
