
// model parameters
var model_params = {
    max_act: 20,
    lambda_act: 200,

    lambda_area: 50,
    lambda_perimeter: 2,
    T: 20,
    J: [[0, null, 20], [null, null, null] ,[20, null, 100]], // index 0='ECM', 1='Skin', 2='Tcell'
    A: [0, 152, 500], // target areas, index 0=ECM, 1=Skin, 2=Tcell
    P: [0, 145, 340], // target perimeneters
}

var visual_params = {
    lattice_width: 800,
    lattice_height: 500,
    n_tissue_cells: 300,
    n_ameboid_cells: 1,

}

var canvasWidth;
var canvasHeight;
const CELL_TYPES = {ECM: 0, SKIN: 1, TCELL: 2};

var lattice_matrix = [[]]
var act_cells = [];

function setup() {
    canvasWidth = 600;
    canvasHeight = 550;
    canvas = createCanvas(canvasWidth, canvasHeight);
    canvas.parent("simulation-canvas");

    // initialize lattice matrix
    for(var i = 0; i < canvasHeight; ++i) {
        lattice_matrix[i] = [];
        for(var j = 0; j < canvasWidth; ++j) {
            lattice_matrix[i][j] = new Site(i, j, null, CELL_TYPES.ECM, 0); // activity value 0 for any medium cell
        }
    }
    // generate n_tissue_cells (skin)
    //cellArea = (canvasWidth * canvasHeight) / n_tissue_cells;


    // seed amoeboid cells
    var tCell1 = seedTCell(1, 1600, 500, 250);
    act_cells.push(tCell1);

    updatePixels();
}

function seedSkinCells(cellArea, canvasWidth, canvasHeight) {
    var b = 0;

}

// draw square Tcell a pos x,y
function seedTCell(cellId, cellArea, x, y) {
    //offset_x = 0;
    //offset_y = 0;
    var cellSide = Math.round(Math.sqrt(cellArea));
    var newCell = new Cell(cellId, cellSide*cellSide, 4*cellSide, null);

    // sites that belong to cell
    for (var i=x; i < x+cellSide; i++) {
        for (var j=y; j < y+cellSide; j++) {
            set(i, j, color(255, 204, 0));
            site = lattice_matrix[i][j];
            site.cellIdentity = cellId;
            site.cellType = CELL_TYPES.TCELL;
            newCell.sites.push(site);
        }
    }
    return newCell;
}

function monteCarloStep() {
    // make canvasWidth*canvasHeight copy attempts
    for (var i=0; i < canvasWidth*canvasHeight; i++) {
        // pick random site u ---- only from broder pixels???
        var u_i = Math.floor((Math.random() * canvasHeight));
        var u_j = Math.floor((Math.random() * canvasWidth));
        var u = lattice_matrix[u_i][u_j];
        // pick random neighbor v of u
        var neighbors = neighbors(u);
        var v = neighbors[Math.floor(Math.random()*neighbors.length)];
        // attempt copy of cell identity only if u and v belong to different cells
        if (u.cellIdentity != v.cellIdentity) { 
            var dH = deltaH_adhesion(u, v) 
                    + deltaH_area(u.cellIdentity, v.cellIdentity) + deltaH_perimeter();
            var dHact = deltaHact(u, v);
            var biased_dH = dH - dHact;
            // probabilistic success
            if (doCopy(biased_dH)) {
                // transformations in site v

            }
        }

    }
}

function hamiltonianGlobal() {
    var adh_term = 0;
    var area_term = 0;
    var perim_term = 0;
    // adhersion term
    for (let cell of act_cells) {
        for (let u of cell.sites) {
            var diffNeighbors = neighborsDiffCellId(u);
            for (let v of diffNeighbors) {
                adh_term += J[u.cellType][v.cellType];
            }
        }
    }
    // area term
    for (let cell of act_cells) {
        var areaDiff = cell.area - model_params.A[cell.sites[0].cellType];
        area_term += (model_params.lambda_area * Math.pow(areaDiff, 2));
    }
    // perimeter term
    for (let cell of act_cells) {
        var perimDiff = cell.perimeter - model_params.P[cell.sites[0].cellType];
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
        adh_term += J[u.cellType][v.cellType];
    }
    return adh_term;
}

// Change in the adhesion term if cell identity of site u is copied to site v.
function deltaH_adhesion(u, v) {
    // shallow copy
    var v_potential = {...v};
    v_potential.cellIdentity = u.cellIdentity;
    v_potential.cellType = u.cellType;
    return adhesionAtSite(v_potential) - adhesionAtSite(v);
}
// Area
function areaConstraintAtCell(cellId, a_change) {
    // Note: medium cells do not have area constraint, so its area contraint term is 0.
    // For the rest, compute it.
    var cell = getCellById(cellId);
    var areaDiff = (cell.area + a_change) - model_params.A[cell.sites[0].cellType];
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
function perimConstraint() {

}

function deltaH_perimeter() {
    return 0;
}

// additional term for the Act model
function deltaHact(u, v) {
    return (model_params.lambda_act / model_params.max_act) * (geomMeanNeighborhood(u) - geomMeanNeighborhood(v));
}

function getCellById(id) {
    var cell = null;
    for (let c of act_cells) {
        if (c.id == id) {
            cell = c; // is this copied by reference???
            break;
        }
    }
    return cell; 
}

function neighbors(site) {
    var neighbors = [];
    if (site.i == 0) { 
        if (site.j == 0) { // top left case
            neighbors.push(lattice_matrix[0][1]);
            neighbors.push(lattice_matrix[1][1]);
            neighbors.push(lattice_matrix[1][0]);
        } else if (site.j == canvasWidth-1) { // top right case
            neighbors.push(lattice_matrix[0][canvasWidth-2]);
            neighbors.push(lattice_matrix[1][canvasWidth-2]);
            neighbors.push(lattice_matrix[1][canvasWidth-1]);
        } else {
            neighbors.push(lattice_matrix[0][site.j-1]);
            neighbors.push(lattice_matrix[0][site.j+1]);
            neighbors.push(lattice_matrix[1][site.j-1]);
            neighbors.push(lattice_matrix[1][site.j+1]);
            neighbors.push(lattice_matrix[1][site.j]);
        }
    } else if (site.i == canvasHeight-1) {
        if (site.j == 0) { // bottom left case
            neighbors.push(lattice_matrix[canvasHeight-2][0]);
            neighbors.push(lattice_matrix[canvasHeight-2][1]);
            neighbors.push(lattice_matrix[canvasHeight-1][1]);
        } else if (site.j == canvasWidth-1) { // bottom right case
            neighbors.push(lattice_matrix[canvasHeight-2][canvasWidth-1]);
            neighbors.push(lattice_matrix[canvasHeight-2][canvasWidth-2]);
            neighbors.push(lattice_matrix[canvasHeight-1][canvasWidth-2]);
        } else {
            neighbors.push(lattice_matrix[canvasHeight-1][site.j-1]);
            neighbors.push(lattice_matrix[canvasHeight-1][site.j+1]);
            neighbors.push(lattice_matrix[canvasHeight-2][site.j-1]);
            neighbors.push(lattice_matrix[canvasHeight-2][site.j+1]);
            neighbors.push(lattice_matrix[canvasHeight-2][site.j]);
        }
    } else if (site.j == 0) {
        neighbors.push(lattice_matrix[site.i-1][0]);
        neighbors.push(lattice_matrix[site.i+1][0]);
        neighbors.push(lattice_matrix[site.i-1][1]);
        neighbors.push(lattice_matrix[site.i+1][1]);
        neighbors.push(lattice_matrix[site.i][1]);
    } else if (site.j == canvasWidth-1) {
        neighbors.push(lattice_matrix[site.i-1][canvasWidth-1]);
        neighbors.push(lattice_matrix[site.i+1][canvasWidth-1]);
        neighbors.push(lattice_matrix[site.i-1][canvasWidth-2]);
        neighbors.push(lattice_matrix[site.i+1][canvasWidth-2]);
        neighbors.push(lattice_matrix[site.i][canvasWidth-2]);
    } else { // internal case
        neighbors.push(lattice_matrix[site.i-1][site.j-1]);
        neighbors.push(lattice_matrix[site.i-1][site.j]);
        neighbors.push(lattice_matrix[site.i-1][site.j+1]);
        neighbors.push(lattice_matrix[site.i][site.j-1]);
        neighbors.push(lattice_matrix[site.i][site.j+1]);
        neighbors.push(lattice_matrix[site.i+1][site.j-1]);
        neighbors.push(lattice_matrix[site.i+1][site.j]);
        neighbors.push(lattice_matrix[site.i+1][site.j+1]);
    }
    return neighbors;
}

function neighborsDiffCellId(site) {
    var diffNeighbors = [];
    var neighbors = neighbors([site.i, site.j]);
    for (let n of neighbors) {
        if (n.cellIdentity != site.cellIdentity) {
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
    var neighbors = neighbors(site);
    var mult = site.activityValue;
    for (let n of neighbors) {
        // only neighbor sites that belong to same cell
        if (site.cellIdentity == n.cellIdentity) {
            mult *= n.activityValue;
        }
    }
    return Math.pow(mult, 1/(neighbors.length+1));
}

function doCopy(dH) {
    if (dH < 0) return true;
    return Math.random() < Math.exp(-dH / model_params.T);
}

function draw() {
    background(235);
    monteCarloStep();
    updatePixels();
}
