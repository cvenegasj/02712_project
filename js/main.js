// main.js

var cpm;
var lattice_width = 200;
var lattice_height = 200;
var dx = 3;
var cellID = 0
var step = 0
var drawBorders = true;


function setup() {
    noStroke();
    canvas = createCanvas(dx*lattice_width, dx*lattice_height);
    canvas.parent("simulation-canvas");
    cpm = new CPM()
    cpm.addCell(CELL_TYPES.TIP, 95, 95, 105, 105)


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
    displayStep()
}


function displayStep() {
    text("Step " + step.toString(), 2, 0.005*dx*lattice_width, 100, 50)
}
