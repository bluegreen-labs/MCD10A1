//----------------------------------------------------------------------------------
//
// Google Earth Engine code to create a combined MODIS snow product, fusing both
// Terra and Aqua MxD10A1 data streams. Both data streams are prone to 
// underestimate snow cover (due to MODIS properties and sky conditions). A fused
// products aleviates these issues partially by taking the maximum value of the
// fractional snow cover bands of each data stream.
//
// This raw data is stored in the maxCollection imageCollection. This data is used
// calculate the date of first snow melt and the date of snow accumulation. These
// dates are determined registering the DOY where fractional snow cover drops below
// or the day before it last exceeds 5%, only earliest and latest instances are noted. 
// A final product can be derived, namely the maximum snow free period, by taking the
// difference between the snow accumulation and snow melt dates.
// 
// Note that this code is writen for the northern hemisphere, and will not work
// for the southern hemisphere without adjustments (sum days between Jan 1. and
// the northern hemisphere snow melt day and from the accumulation day until 
// Dec. 31th to get to the snow free period, and obviously reverse the captions
// on the dates melt becomes accumulation and vice versa).
//
// Code written by Koen Hufkens, April 2016.
// 
//----------------------------------------------------------------------------------

// -- settings --

// years to process (from start year t0 to end year t1)
// [2003 is the first year with full coverage]
var t0 = "2003", t1 = "2015";

// -- functions --

var maxVal = function(feature) {
  // convert feature to true image
  // for some reason the join messes this up
  var img = ee.Image(feature);
  var b1 = img.select('Fractional_Snow_Cover');
  var b1m = b1.updateMask(b1.lt(100));
  var b2 = img.select('Fractional_Snow_Cover_1');
  var b2m = b2.updateMask(b2.lt(100));
  var startTime = b1.get('system:time_start');
  // take maximum value of extracted bands, retain the
  // start time, but unmask the masked values as not to
  // mess with the DOY routine.
  var maxval = b1m.max(b2m).set('system:time_start',startTime).unmask(200);
  return maxval;
};

var difference = function(feature) {
  var img = ee.Image(feature);
  var b1 = img.select('snowmelt');
  var b2 = img.select('snowacc');
  var startTime = b1.get('system:time_start');
  return b2.subtract(b1).set('system:time_start',startTime);
};

// set the snow melt date
var returnDoy = function(feature) {
  var image = ee.Image(feature);
  // Compute time of the image in fractional years relative to the start.
  var doy = ee.Date(image.get('system:time_start')).difference(start, 'day').add(1);
  // replace values smaller than 5% snow cover with
  // their respective DOY
  var replacement = image.expression(
      'SNOW <= 5 ? DOY : 9999', {
        'SNOW': image.select('Fractional_Snow_Cover'),
        'DOY': doy
  }).toInt();
  
  // mask anything above 5% coverage as empty
  var withObs = replacement.neq(9999);
  return replacement.updateMask(withObs);
};

// stick two bands together which are linked with an
// inner join
var concatBands = function(feature) {
  return ee.Image.cat(feature.get('primary'),feature.get('secondary'));
};

// This function adds a time band to the image,
// for subsequent regression analysis
var createTimeBand = function(feature) {
  var image = ee.Image(feature);
  return image.addBands(image.metadata('system:time_start'));
};

// Specify an equals filter for image timestamps.
var filterTimeEq = ee.Filter.equals({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

// -- main -- 

// If a longer time series is selected, calculate statistics
// for each year independently. This can probably be done
// with some sort of reduce statement, but the documentation is a bit
// of a mess when it comes to custom functions.
for (var i = parseInt(t0); i <= parseInt(t1); ++i) {
  
  // int to string
  var evalYear = i.toString();
  
  // Load MODIS Fractional Snow Cover product (aqua / terra)
  var terraCollection = ee.ImageCollection('MOD10A1')
    .select('Fractional_Snow_Cover')
    .filterDate(evalYear.concat("-01-01"),evalYear.concat("-12-31"))
    .sort('system:time_start', true);
    //.filterBounds(NorthSlope);
  
  var aquaCollection = ee.ImageCollection('MYD10A1')
    .select('Fractional_Snow_Cover')
    .filterDate(evalYear.concat("-01-01"),evalYear.concat("-12-31"))
    .sort('system:time_start', true);
    //.filterBounds(NorthSlope);
  
  // create a collection of both MOD/MYD data matched
  // joined per band using an inner join
  var innerJoin = ee.Join.inner();
  
  // Apply the join.
  var innerJoinedCollection = innerJoin.apply(terraCollection, aquaCollection, filterTimeEq);
  
  // combine matching images into one image with (2) bands
  var joinedCollection = innerJoinedCollection.map(concatBands);
  
  // calculate the maximum value across the two bands MOD/MYD
  var yearlyCollection = joinedCollection.map(maxVal);

  // how many bands in the dataset?
  // should be close to 365
  // print(yearlyCollection.size());

  // Get the starting time reference.
  var start = ee.Date(yearlyCollection.first().get('system:time_start'));

  // subsitute all values < 5% with their respective DOY
  var doyCollection = yearlyCollection.map(returnDoy);

  // calculate the snow melt and snow fall (accumulation)
  var snowMeltAcc = ee.ImageCollection(doyCollection)
    .reduce(ee.Reducer.minMax()).unmask(366);
    
  // The latter part is rather ugly, not optimal (but works)
  // if first in loop swap with the original layer
  if (evalYear == t0){

    // snow melt
   var tmpImg = snowMeltAcc
    .select(['constant_min'],['snowmelt'])
    .set('system:time_start',i);
   var snowMelt = ee.ImageCollection(tmpImg);
   
   // snow accumulation
   var tmpImg = snowMeltAcc
    .select(['constant_max'],['snowacc'])
    .set('system:time_start',i);
   var snowAcc = ee.ImageCollection(tmpImg);
   
  }else{
    // add snow melt layers
    var tmpImg = snowMeltAcc.select(['constant_min'],['snowmelt'])
      .set('system:time_start',i); 
    var snowMelt = snowMelt.merge(ee.ImageCollection(tmpImg));
    
  // add snow accumulation layers
    var tmpImg = snowMeltAcc.select(['constant_max'],['snowacc'])
      .set('system:time_start',i); 
    var snowAcc = snowAcc.merge(ee.ImageCollection(tmpImg));
  }
}

// calculate the snow free season length
// in number of days (for masking or further analysis)

// Apply the join
var innerJoinedCollection = innerJoin.apply(snowMelt, snowAcc, filterTimeEq);

// combine matching images into one image with (2) bands
var joinedCollection = innerJoinedCollection.map(concatBands);

// take the difference of the two bands
var snowFree = joinedCollection.map(difference);

// calculate median snow free season length and
// and create a mask from the result
var snowFreeMedian = ee.ImageCollection(snowFree)
  .median();
  
var snowFreeMask = ee.Image(snowFreeMedian)
  .lt(306);

// -- quick regression analysis -- 

// add a time band
var snowMelt = snowMelt.map(createTimeBand);

// linear regression on the data
var linearFit = ee.ImageCollection(snowMelt).select(['system:time_start', 'snowmelt'])
  .reduce(ee.Reducer.linearFit());

// -- plotting results --

// land-water mask for visualization, removes faulty cover
// values for the oceans and water bodies
var hansenImage = ee.Image('UMD/hansen/global_forest_change_2013');
var data = hansenImage.select('datamask');
var waterMask = data.eq(1);

// final mask, excluding water bodies and
// areas with a too short of a winter season
var mask = waterMask.multiply(snowFreeMask);

// Plot the snow free trend
Map.addLayer(linearFit
  .mask(mask),
  {min: 0, max: [-2, 0, 2], bands: ['scale', 'offset', 'scale']}, 'fit');