
{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "box",
size: [40000.0, 40000.0, 40000.0], // mm, half-height
material: "rock",
color: [0.67, 0.29, 0.0],
invisible: 1,
}

{
name: "GEO",
index: "outerWP",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world",
type: "tube",
r_max: 21750, // radius, mm
size_z: 22000, // half-height, mm
material: "water",
color: [1.0, 1.0, 1.0],
invisible: 1,
}

// Do not put the SS supppor for now
//{
//name: "GEO",
//index: "pmt_support",
//valid_begin: [0, 0],
//valid_end: [0, 0],
//mother: "outerWP",
//type: "sphere",
//r_min: 19930.0, // Inner radius, mm 
//r_max: 28180.0 // outer radius, mm
///////// Need to check the SS size
//material: "stainless_steel",
//color: [1.0, 1.0, 1.0],
//invisible: 1,
//}


{
name: "GEO",
index: "innerPMT",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "outerWP",
type: "pmtarray",
pmt_model: "r1408", 
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner",
efficiency_correction: 1.000,  
add_concentrator: 0, // Flag: 0 = no concentrators, 1 = concentrators
pos_table: "PMTINFO_CD_LPMT_JUNO",
orientation: "point", // Aim all PMTs at a point
	     	      // "manual" means there is a dir_x, dir_y, dir_z
	              // in pos_table for manual orientation of PMTs
orient_point: [0.0, 0.0, 0.0],
}


{
name: "GEO",
index: "WPPMT",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "outerWP",
type: "pmtarray", 
pmt_model: "r1408",
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner",
efficiency_correction: 1.000,
add_concentrator: 0,
pos_table: "PMTINFO_WP_LPMT_JUNO",
orientation: "manual",
invisible: 1,
}


{
name: "GEO",
index: "WPWALPMT",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "outerWP",
type: "pmtarray", 
pmt_model: "r1408",
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner",
efficiency_correction: 1.000,
add_concentrator: 0,
pos_table: "PMTINFO_WP_WALPMT_JUNO",
orientation: "manual",
}

{
name: "GEO",
index: "acrylic",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "outerWP",
type: "sphere",
r_min: 17700, // Inner radius, mm
r_max: 17824, // Outer radius, mm
material: "acrylic_sno",
color: [0.0, 0.3, 1.0, 0.2],
}

{
name: "GEO",
index: "LS",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "outerWP",
type: "sphere",
r_max: 17700, 
material: "scintillator",
color: [0.2, 0.1, 0.1, 0.5],
}
