const CELL_TYPES = {
  ECM: 0,
  SKIN: 1,
  TIP: 2,
  TUMOR:3,
  HYPOXIC: 4
}

const COLOR = [
  "#EEEEEE",
  "#FF6E6E",
  "#FFCC00",
  "#77BBFF",
  "#4C76B5"
]

// model parameters
var model_params = {
  max_act: 40,
  lambda_act: 200,

  lambda_area: 50,
  lambda_perimeter: 2,
  T: 20,
  J: [[0, null, 20], [null, null, null] ,[20, null, 100]], // index 0='ECM', 1='Skin', 2='Tip cell'
  A: [0, 152, 500], // target areas, index 0=ECM, 1=Skin, 2=Tip cell
  P: [0, 145, 340], // target perimeneters
}

var now_active = [];

class CPM {
  constructor() {
    this.cells = []
    this.lattice_matrix = [[]]
    this.active_sites = []
    for(var i = 0; i < lattice_height; ++i) {
      this.lattice_matrix[i] = [];
      for(var j = 0; j < lattice_width; ++j) {
          this.lattice_matrix[i][j] = new Site(i, j, null, CELL_TYPES.ECM, 0); // activity value 0 for any medium cell
      }
    }
  }

  initialize() {
    
  }

  addCell(cellType, x0, y0, x1, y1) {
    var sites = [], borderSites = [];
    for (var i=x0; i < x1; i++) {
        for (var j=y0; j < y1; j++) {
            var site = this.lattice_matrix[i][j];
            site.cellId = cellID;
            site.cellType = cellType;
            sites.push(site);
            if (i==x0 || i==x1-1 || j==y0 || j==y1-1) {
              borderSites.push(site);
            }
        }
    }
    this.cells.push(new Cell(cellType, cellID, sites, borderSites));
    cellID++;
  }
  
  monteCarloStep() {
    // print(cpm.cells[0].sites)
    // print(cpm.cells[0].area())
    for (var i=0; i < lattice_width*lattice_height; i++) {
      // pick random site u ---- only from broder pixels???
      var u_i = Math.floor((Math.random() * lattice_height));
      var u_j = Math.floor((Math.random() * lattice_width));
      var u = this.lattice_matrix[u_i][u_j];
      // pick random neighbor v of u
      var neighbors = neighborsOf(u);
      var v = neighbors[Math.floor(Math.random()*neighbors.length)];
      // attempt copy of cell identity only if u and v belong to different cells
      if (u.cellId != v.cellId) { 
        var dH = deltaH_adhesion(u, v) 
        + deltaH_area(u.cellId, v.cellId) + deltaH_perimeter();
        var dHact = deltaHact(u, v);
        var biased_dH = dH - dHact;
        // probabilistic success
            if (doCopy(biased_dH)) {
                // print(u)
                // print(v)
                // print(cpm.cells[0].sites)
                siteUpdate(u, v)// transformations in site v
            }
        }
    }
    this.updateActiveSites();
    step++;
  }

  updateActiveSites() {
    for (let i=0; i<this.active_sites.length; i++) {
      if (this.active_sites.activityValue == 0) {
        this.active_sites.splice(i, 1);
      } else {
        this.active_sites.activityValue--;
      }
    }
    for (let v of this.active_sites) {
      v.activityValue=model_params.max_act;
    }
    this.active_sites.push(now_active);
    now_active = [];
  }
}


function hamiltonianGlobal() {
  var adh_term = 0;
  var area_term = 0;
  var perim_term = 0;
  // adhersion term
  for (let cell of cpm.cells) {
      for (let u of cell.sites) {
          var diffNeighbors = neighborsDiffCellId(u);
          for (let v of diffNeighbors) {
              adh_term += J[u.cellType][v.cellType];
          }
      }
  }
  // area term
  for (let cell of cpm.cells) {
      var areaDiff = cell.area() - model_params.A[cell.sites[0].cellType];
      area_term += (model_params.lambda_area * Math.pow(areaDiff, 2));
  }
  // perimeter term
  for (let cell of cpm.cells) {
      var perimDiff = cell.perimeter() - model_params.P[cell.sites[0].cellType];
      perim_term += (model_params.lambda_perimeter * Math.pow(perimDiff, 2));
  }
  var H = adh_term + area_term + perim_term;
  return H;
}

/* Hamiltonian terms computed at single lattice sites */
// Adhesion 
function adhesionAtSite(u) {
  var adh_term = 0;
  var diffNeighbors = neighborsDiffCellId(u);
  for (let v of diffNeighbors) {
      adh_term += model_params.J[u.cellType][v.cellType];
  }
  return adh_term;
}

// Change in the adhesion term if cell identity of site u is copied to site v.
function deltaH_adhesion(u, v) {
  // shallow copy
  var v_potential = {...v};
  v_potential.cellId = u.cellId;
  v_potential.cellType = u.cellType;
  return adhesionAtSite(v_potential) - adhesionAtSite(v);
}
// Area
function areaConstraintAtCell(cellId, a_change) {
  // Note: medium cells do not have area constraint, so its area contraint term is 0.
  // For the rest, compute it.
  if (cellId == null) {
      return 0;
  }
  var cell = getCellById(cellId);
  var areaDiff = (cell.area() + a_change) - cell.targetArea;
  return model_params.lambda_area * areaDiff * areaDiff;
}

function deltaH_area(srcCellId, tgtCellId) {
  // change due to volume increase in source cell
  var deltaH = areaConstraintAtCell(srcCellId, 1) - areaConstraintAtCell(srcCellId, 0);
  // change due to volume decrease in target cell
  deltaH += areaConstraintAtCell(tgtCellId, -1) - areaConstraintAtCell(tgtCellId, 0);
  return deltaH;
}

// Perimeter
function perimeterConstraintAtCell(cellId, p_change) {
  var cell = getCellById(cellId);
  if (cellId == null) {
    return 0;
  }
  var perimeterDiff = (cell.perimeter() + p_change) - model_params.P[cell.type];
  return model_params.lambda_perimeter * perimeterDiff * perimeterDiff
}

function deltaH_perimeter(srcCellId, tgtCellId) {
  // change in source cell
  var deltaH = perimeterConstraintAtCell(srcCellId, 1) - perimeterConstraintAtCell(srcCellId, 0);
  // change in target cell
  deltaH += perimeterConstraintAtCell(tgtCellId, -1) - perimeterConstraintAtCell(tgtCellId, 0);
  return deltaH;
}

// additional term for the Act model
function deltaHact(u, v) {
  return (model_params.lambda_act / model_params.max_act) * (geomMeanNeighborhood(u) - geomMeanNeighborhood(v));
}

function getCellById(id) {
  for (let c of cpm.cells) {
      if (c.id == id) {
          return c;
      }
  }
}

function neighborsOf(site) {
  var neighbors = [];
  if (site.i == 0) { 
      if (site.j == 0) { // top left case
          neighbors.push(cpm.lattice_matrix[0][1]);
          neighbors.push(cpm.lattice_matrix[1][1]);
          neighbors.push(cpm.lattice_matrix[1][0]);
      } else if (site.j == lattice_width-1) { // top right case
          neighbors.push(cpm.lattice_matrix[0][lattice_width-2]);
          neighbors.push(cpm.lattice_matrix[1][lattice_width-2]);
          neighbors.push(cpm.lattice_matrix[1][lattice_width-1]);
      } else {
          neighbors.push(cpm.lattice_matrix[0][site.j-1]);
          neighbors.push(cpm.lattice_matrix[0][site.j+1]);
          neighbors.push(cpm.lattice_matrix[1][site.j-1]);
          neighbors.push(cpm.lattice_matrix[1][site.j+1]);
          neighbors.push(cpm.lattice_matrix[1][site.j]);
      }
  } else if (site.i == lattice_height-1) {
      if (site.j == 0) { // bottom left case
          neighbors.push(cpm.lattice_matrix[lattice_height-2][0]);
          neighbors.push(cpm.lattice_matrix[lattice_height-2][1]);
          neighbors.push(cpm.lattice_matrix[lattice_height-1][1]);
      } else if (site.j == lattice_width-1) { // bottom right case
          neighbors.push(cpm.lattice_matrix[lattice_height-2][lattice_width-1]);
          neighbors.push(cpm.lattice_matrix[lattice_height-2][lattice_width-2]);
          neighbors.push(cpm.lattice_matrix[lattice_height-1][lattice_width-2]);
      } else {
          neighbors.push(cpm.lattice_matrix[lattice_height-1][site.j-1]);
          neighbors.push(cpm.lattice_matrix[lattice_height-1][site.j+1]);
          neighbors.push(cpm.lattice_matrix[lattice_height-2][site.j-1]);
          neighbors.push(cpm.lattice_matrix[lattice_height-2][site.j+1]);
          neighbors.push(cpm.lattice_matrix[lattice_height-2][site.j]);
      }
  } else if (site.j == 0) {
      neighbors.push(cpm.lattice_matrix[site.i-1][0]);
      neighbors.push(cpm.lattice_matrix[site.i+1][0]);
      neighbors.push(cpm.lattice_matrix[site.i-1][1]);
      neighbors.push(cpm.lattice_matrix[site.i+1][1]);
      neighbors.push(cpm.lattice_matrix[site.i][1]);
  } else if (site.j == lattice_width-1) {
      neighbors.push(cpm.lattice_matrix[site.i-1][lattice_width-1]);
      neighbors.push(cpm.lattice_matrix[site.i+1][lattice_width-1]);
      neighbors.push(cpm.lattice_matrix[site.i-1][lattice_width-2]);
      neighbors.push(cpm.lattice_matrix[site.i+1][lattice_width-2]);
      neighbors.push(cpm.lattice_matrix[site.i][lattice_width-2]);
  } else { // internal case
      neighbors.push(cpm.lattice_matrix[site.i-1][site.j-1]);
      neighbors.push(cpm.lattice_matrix[site.i-1][site.j]);
      neighbors.push(cpm.lattice_matrix[site.i-1][site.j+1]);
      neighbors.push(cpm.lattice_matrix[site.i][site.j-1]);
      neighbors.push(cpm.lattice_matrix[site.i][site.j+1]);
      neighbors.push(cpm.lattice_matrix[site.i+1][site.j-1]);
      neighbors.push(cpm.lattice_matrix[site.i+1][site.j]);
      neighbors.push(cpm.lattice_matrix[site.i+1][site.j+1]);
  }
  return neighbors;
}

function neighborsDiffCellId(site) {
  var diffNeighbors = [];
  var neighbors = neighborsOf(site);
  for (let n of neighbors) {
      if (n.cellId != site.cellId) {
          diffNeighbors.push(n);
      }
  }
  return diffNeighbors;
}

function geomMeanNeighborhood(site) {
  // gemoetric mean of any site that belongs to medium is 0
  if (site.cellType == CELL_TYPES.ECM) {
      return 0;
  }
  var neighbors = neighborsOf(site);
  var mult = site.activityValue;
  for (let n of neighbors) {
      // only neighbor sites that belong to same cell
      if (site.cellId == n.cellId) {
          mult *= n.activityValue;
      }
  }
  return Math.pow(mult, 1/(neighbors.length+1));
}

function doCopy(dH) {
  if (dH < 0) return true;
  return Math.random() < Math.exp(-dH / model_params.T);
}

// state of site v gets updated to state of site u 
function siteUpdate(u, v) {
  var cell1 = getCellById(u.cellId); 
  var cell2 = getCellById(v.cellId);
  // v.cellType = CELL_TYPES.SKIN;
  now_active.push(v);
  if (cell1 != null) {
    cell1.addSite(v);
  }
  else {
    v.cellId = null;
    v.cellType = 0;
  }
  if (cell2 != null) {
    cell2.removeSite(v);
  }
}