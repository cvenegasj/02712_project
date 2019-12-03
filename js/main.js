// main.js

var cpm;
var lattice_width = 200;
var lattice_height = 200;
var dx = 3;
var dt = 1;
var cellID = 0
var step = 0
var drawBorders = true;
var updateConcentrations = false;


function setup() {
    noStroke();
    canvas = createCanvas(dx*lattice_width, dx*lattice_height);
    canvas.parent("simulation-canvas");
    cpm = new CPM()
    // cpm.addCell(CELL_TYPES.TIP, 95, 95, 105, 105)
    // cpm.addCell(CELL_TYPES.TIP, 75, 75, 85, 85)
    // cpm.initialize(SETTING.CHEMOTAXIS)
    cpm.initialize(SETTING.TUMOR)

    // background(235);
    // for (let c of cpm.cells) {
    //     c.draw()
    // }
    // cpm.monteCarloStep();
}


function draw() {
    background(235);
    for (let c of cpm.cells) {
        c.draw()
    }
    cpm.monteCarloStep();
    if (updateConcentrations) {
        cpm.updateConcentrations();
    }
    cpm.growth()
    displayStep()
}


function displayStep() {
    text("Step " + step.toString(), 2, 0.005*dx*lattice_width, 100, 50)
}
