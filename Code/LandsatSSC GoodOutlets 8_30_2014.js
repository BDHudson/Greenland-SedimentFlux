// LandsatSSC GoodOutlets 8_30_2014
// ----------------------------------------------------------
// Ben Hudson, BDHudson@Gmail.com
// started 8-19-2014
// PURPOSE:
// To look at a number of rivers in Greenland, calculate the Suspended 
// Sediment Concentration in them, and the number of times an area has been 
// inundated (water was present) and extract it to a CSV file. 
// ----------------------------------------------------------
// ----------------------------------------------------------
// START FUNCTIONS 
// ----------------------------------------------------------
// ----------------------------------------------------------
// FUNCTION 1 - This function makes a ROI box and returns 
function makeBox(ROI_lat, ROI_lon,dn,de) {

  //Position, decimal degrees
  //Earth’s radius, sphere
  var EarthRadius = 6378137;
  
  //Coordinate offsets in radians
  var dLat = dn/EarthRadius;
  var dLon = de/(EarthRadius*Math.cos(Math.PI*ROI_lat/180));
  
  //OffsetPosition, decimal degrees
  var latPlus = ROI_lat + dLat * 180/Math.PI;
  var lonPlus = ROI_lon + dLon * 180/Math.PI;
  
  //OffsetPosition, decimal degrees
  var latMinus = ROI_lat - dLat * 80/Math.PI;
  var lonMinus = ROI_lon - dLon * 80/Math.PI;
  
  //print('dn:' + dn);
  //print('de:' + de);

  //print('dLat:' + dLat);
  //print('dLon:' + dLon);
  //print('latPlus:' + latPlus);
  //print('lonPlus:' + lonPlus);
  //print('latMinus:' + latMinus);
  //print('lonMinus:' + lonMinus);
  
  //ALWAYS COUNTER CLOCKWISE!!!
  // Create a new feature collection by defining geometries and 
  // associated properties.
  var fc = new ee.FeatureCollection([
      ee.Feature(
          ee.Geometry.Polygon(
              [[lonPlus,latPlus],[lonMinus,latPlus],[lonMinus,latMinus],[lonPlus,latMinus]]),
          {'name': 'ROI', fill: 1})
      
      ]);
  return fc;    
  }
  
// ----------------------------------------------------------
// FUNCTION 2 - PROCESS SSC  

function processSSC(image) {
  // Calculate SSC 
  var SSC = image.expression('4.9 * exp(38.0 * band4)',
      { band4: image.select('B4')});
  // Create Mask 1
  var water = image.expression('band5 < .05 ? 1 : 0 ',
      { band5: image.select('B5')
      });
  // Create mask 2 - this is a snow/ice mask
  var band2mask = image.expression('band2 > .3 ? 0 : 1 ',
      { band2: image.select('B2')
      });
  // Create mask 3 - this is a shadow mask    
  var band1mask = image.expression('band1 > .15 ? 1 : 0 ',
      { band1: image.select('B1')
      });
  // Create mask 4 - this is a ice 
  var band4mask = image.expression('band4 <= .25 ? 1 : 0 ',
      { band4: image.select('B4')
      });
  // Create mask 5 - this is a ice 
  var band3mask = image.expression('band3 <= .28 ? 1 : 0 ',
      { band3: image.select('B3')
      });

  // Combine masks 
  var allMasks_1 = band1mask.multiply(water).float();
  var allMasks_2 = allMasks_1.multiply(band2mask).float();
  var allMasks_3 = allMasks_2.multiply(band4mask).float();
  var allMasks = allMasks_3.multiply(band3mask).float();
 
  var combo = SSC.multiply(allMasks);
  
  var justSSC = combo.mask(allMasks);
  allMasks = allMasks.mask(allMasks);
  
  // This returns the SSC and the innundation maps
  image = image.addBands([allMasks,justSSC]);
  //image = image.addBands([justSSC,allMasks]);
  
  // Return the updated image.
  return image;
   
  }

// ----------------------------------------------------------
// FUNCTION 3 - mask out entire collection
// https://sites.google.com/site/earthengineapidocs/javascript-api/javascript-developers-guide/image-methods
function maskCollection(collection,maskImage) {
  
  // text 
  var cMasked = collection.map(function(image) {
  
  // throw out pixels in the sum image that were not greater than or equal to x 
  
  
  //maskedImage = image.maskImage.gte(10);
  var mi = ee.Image(image).mask(maskImage.gte(10));
  //print('masked Image:', maskImage.getInfo());
  
  return mi;  
    
  });

  return cMasked; 
  
  } 


// ----------------------------------------------------------
// FUNCTION 4 - create CSV 

function roiCSV(extractCollection,roi,on) {

  // Select only the band you want to extract - for simplicity leaving it to one band at a time. 

  var reducer = ee.Reducer.mean()
              .combine(ee.Reducer.min(), '', true)
              .combine(ee.Reducer.max(), '', true)
              .combine(ee.Reducer.median(), '', true)
              .combine(ee.Reducer.sum(),'',true);
              
  var c4 = ee.FeatureCollection(extractCollection.map(function(image) {
  
    var reduced = image.reduceRegion({
      reducer: reducer,
      geometry: roi,
      crs: 'EPSG:32622',
      scale: 30
    });
    // null below means no geometry is attached 
    var result = ee.Feature(null, reduced);
    
    var newResult = result.set('fName', image.get('LANDSAT_SCENE_ID'));
    
    return newResult.set('date', image.get('DATE_ACQUIRED'));
    
  }));
  
  // exclude scenes that had no valid pixels in the roi
  // then sorts by date 
  
  c4 = c4.filter(ee.Filter.neq('constant_min', null)).sort('date');
  
    // Add index number 
    var addindex = function(feature) {
  
    var newFeature = feature.set({'indexNo': i});
    return newFeature;
    };

    var c5 = c4.map(addindex);
   
    // Add Outlet number 
    var addOutlet = function(feature) {
  
    var newFeature = feature.set({'outletNo': on});
    return newFeature;
    };

    var c6 = c5.map(addOutlet); 
    
    // Add File Name / scene id
    //var addFN = function(feature) {
  
    //var newFeature = feature.set('fName', image.get('LANDSAT_SCENE_ID'));
    //return newFeature;
    //};

    //var c7 = c6.map(addFN); 
    
  return c6 
  
  //print(c4.getInfo());
  //print(roi.getInfo());
  
  //print(c4.getDownloadURL('CSV','outletNo,indexNo,date,constant_mean,constant_min,constant_max,constant_median,constant_sum'));
  } // end 3b

// ----------------------------------------------------------
// FUCTION - add 

// ----------------------------------------------------------
// START MAIN BODY OF CODE
// ----------------------------------------------------------

// ----------------------------------------------------------
// INPUT - FOR ROIs

var ROI_lat = [68.20777,68.13102,67.701011,67.51,67.185415,
              67.061317,67.082267,66.976691,66.963,65.1835,
              64.450545,64.087475,63.534469,63.4300,65.689969,
              65.730,65.79821,65.560,66.231457,66.542818,
              66.585785,66.636804,66.7040,66.781086,67.297417, // 20 24
              67.435,67.595,67.661205,67.781153,67.922561, // 25 29
              67.882076,67.871461,68.060591,72.340,72.670,// 30 34
              63.2855,71.458313,72.070,72.564436,74.314,// 35 39
              74.242419,74.160,75.806797,77.78997,73.910052,// 40 44
              74.119208,79.604196,79.665,79.950273,80.19440,// 45 49
              80.44803,80.552,80.600679,80.576788,80.747081];// 50 55

var ROI_lon = [-50.595339,-50.48528,-49.977037,-50.056,-50.384126,
              -50.198005,-50.253737,-50.067177,-50.003568,-50.70,
              -49.681245,-49.985373,-50.721176,-49.698964,-50.569981,
              -50.712657,-51.067508,-50.097295,-49.847442,-49.642927,
              -49.617716,-49.578545,-49.867,-49.834995,-49.98, // 20 24
              -49.815,-49.901784,-49.871,-49.904381,-50.2,// 25 29
              -50.104137,-50.064264,-50.350563,-53.697092,-54.065,// 30 34
              -49.575,-28.875,-26.431631,-26.510087,-22.481764,// 35 39
              -22.380453,-22.345,-22.349413,-21.69594,-24.106265,// 40 44
              -24.467952,-23.27684,-23.648715,-23.82,-24.01749,// 45 49
              -25.441823,-25.656082,-26.121851,-26.737433,-26.402776];// 50 55

//offsets in meters
var dn = [1500,1500,1000,1500,1700, // 0 4
          500,500,500,600,1500, // 5 9
          1000,1500,3000,2000,1500, // 10 14
          3000,1500,1250,500,1000, // 15 19
          1100,200,1500,1000,1500,// 20 24
          1250,1000,500,1200,1300,// 25 29
          600,500,500,2200,750,
          750,1000,1500,750,200,// 35 39
          500,500,500,150,350,// 40 44
          1000,750,150,650,400,// 45 49
          1250,200,750,200,1500];// 50 55


var de = [1500,1500,1500,500,1600, // 0 4
          900,1000,500,500,1500, // 5 9
          1000,1500,1500,250,1500,
          1500,1500,250,500,1250,
          500,200,1000,1000,1500,// 20 24
          800,650,750,1200,2000,// 25 29
          600,500,500,1500,750,// 30 34
          1000,450,1000,750,200,// 35 39
          500,500,750,150,350,// 40 44
          1000,500,150,650,400,// 45 49
          1250,200,750,600,1500];// 50 55

var outletNumber = [2,3,5,6,7,
                    8,9,10,11,12,
                    14,15,16,17,24,
                    25,26,28,30,31,
                    32,33,34,35,37,
                    38,40,42,43,45,
                    46,47,48,54,55,
                    62,102,103,106,109,
                    110,112,117,125,129,
                    131,137,139,140,143,
                    155,157,158,159,160];

// INPUT - FOR ROIs

//var ROI_lat = [68.342164,68.3145,68.28929,68.27401,
//              67.998472,67.93329,67.9210,67.922224,62.622759,
//62.629664,62.485326,62.446853,76.739226,76.717452,
//              76.700007];
//,
//var ROI_lon = [-51.182283,-51.212,-51.18836,-51.132799,
//              -50.385951,-50.245839,-50.303719,-50.28595,-50.106561,
//              -50.170143,-50.232358,-49.98939,-23.041698,-23.15,
//              -23.080];

//offsets in meters
//var dn = [300,300,250,500,
//          500,400,200,200,250,
//          750,1500,1500,750,500,
//          1000];
//var de = [300,300,250,500,
//          500,500,200,200,250,
//          500,1500,1500,750,1000,
//          300];// 50 55

//var outletNumber = [1.1,1.2,1.3,1.4,
//                    4.1,4.2,4.3,4.4,18.1,
//                    18.2,18.3,18.4,119.1,119.2,
//                    119.3];
// ----------------------------------------------------------
// INPUT - FOR IMAGE PROCESSING
var dateStart = '1999-06-15';
var dateEnd = '2014-8-15';


// PRE allocate collection that will get merged on each loop 

var BigC4 = new ee.FeatureCollection([
    ee.Feature(null)
    ]);
// ----------------------------------------------------------
// START FOR LOOP THAT DOES MOST OF THE WORK 
// ----------------------------------------------------------
var ii;
//var i =0;

var iRange = [50,51,52,53,54];
for (ii in iRange) {
//for (i in ROI_lat) {
    
    print(ii);
    // js doesn't know i want the value, not just the index
    //this reads the value from the index 
    var i = iRange[ii];
    print(i);
    // ROI MADE 
    var fc = makeBox(ROI_lat[i], ROI_lon[i],dn[i],de[i]);
    //centerMap(-50.67637,67.00758, 10);
    //centerMap(ROI_lon[i],ROI_lat[i], 10);
    
    // COLLECTION PROCESSSED WITHIN THE TINY ROI
    
    var collection = ee.ImageCollection('LE7_L1T_TOA')
       .filterDate(dateStart, dateEnd)
       .filterBounds(fc)
       .filter(ee.Filter.lt('CLOUD_COVER', 50.0))
       .filter(ee.Filter.gt('google:registration_count', 1000)) // These calls check the reistration
       .filter(ee.Filter.gt('google:registration_offset_x', -100))
       .filter(ee.Filter.lt('google:registration_offset_x', 100))
       .filter(ee.Filter.lt('google:registration_offset_y', 100))
       .filter(ee.Filter.gt('google:registration_offset_y', -100));
       
    //var doyFilter = ee.Filter.dayOfYear(190, 210)
    var doyFilter = ee.Filter.dayOfYear(160, 240);
    
    // Apply the filter using the filter method.
    var collection = collection.filter(doyFilter);
    
    
    // Now create a new collection which contains all the same images as the
    // old one, but where each additional bands 
    
    var newCollection = collection.map(processSSC);
    
    // IMAGE LEVEL WORK
    
    // Calculates the median for the image 
    var medianImage = newCollection.median();
    var sumImage = newCollection.sum();

    // Make it into simple, 1 band images 
    // select constant (SSC map)
    // cast it as a unit16 ( 65,535 values)
    // mask where the water mask doesn't have at least 10 
    
    //SSC
    //medianImage = ee.Image(medianImage).select('constant').toUint16().mask(sumImage.select('constant_1').gte(10));
    //sumImage = ee.Image(sumImage).select('constant_1').toUint8().mask(sumImage.select('constant_1').gte(10));

    // mask where the water mask doesn't have at least 10 
    // INUN
    medianImage = ee.Image(medianImage).select('constant_1').toUint16().mask(sumImage.select('constant').gte(10));
    sumImage = ee.Image(sumImage).select('constant').toUint8().mask(sumImage.select('constant').gte(10));
    
    // constant is SSC, constant_1 is innundation
    var extractCollection = newCollection.select('constant');

    // This step gets rid of pixels in each Roi that did not have at least 10 images contributing
    var maskedNewCollection = maskCollection(extractCollection,sumImage); 
    
    //  extract CSV
    var c4 = roiCSV(maskedNewCollection,fc,outletNumber[i]);
    
    // merge c4 and bigC4 in attempt to concatenate all the results 
    BigC4 = BigC4.merge(c4); 
    
    //addToMap(fc);
    //addToMap(medianImage)
  }

//print(c4.getInfo());
print(BigC4.getDownloadURL('CSV','outletNo,indexNo,fName,date,constant_mean,constant_min,constant_max,constant_median,constant_sum','GoodOutlets_8_30_2014_INUN'));

//print(BigC4.getDownloadURL('CSV','outletNo,indexNo,date,constant_mean,constant_min,constant_max,constant_median,constant_sum','GoodOutlets_8_29_2014_INUN'));


// ----------------------------------------------------------
// END FOR LOOP THAT DOES MOST OF THE WORK
// ----------------------------------------------------------
//print(collection.getInfo());  

// PUT ON GOOGLE MAP 

//var SSC_PALETTE = '0000ff, 00ffff,ffff00,ffA500,ff0000,800000';
//var trend_PALETTE = '0000ff, 00ffff,ffff00,ffA500,ff0000,800000'; //'0000ff, ffffff,ff0000';
//var sum_PALETTE = 'ffffff, 000080'; //'0000ff, ffffff,ff0000';

//(medianImage, {'min': 0, 'max':20000,'palette': SSC_PALETTE});
//
//addToMap(sumImage, {'min': 0, 'max':150,'palette': sum_PALETTE})
