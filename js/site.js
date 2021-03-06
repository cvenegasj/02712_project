class Site {

    constructor(i, j, cellId, cellType, activityValue) {
        this.i = i;
        this.j = j;
        this.cellId = cellId;
        this.cellType = cellType; // integer
        this.activityValue = activityValue;
        this.VEGF = 0;
        if (cellType == CELL_TYPES.SKIN) {
            this.O2 = model_params.max_O2;
        } else {
            this.O2 = 0;
        }
    }

    isOnBorder(cellId) {
        var neighbors = neighborsOf(this)
        for (let v of neighbors) {
            if (v.cellId != this.cellId) {
                return true
            }
        }
        return false
    }

    draw() {
        push();
        // fill(color(COLOR[this.cellType]));
        // var a = this.activityValue/model_params.max_act;
        // if (a > 0.5) {
            //     fill(Math.floor(256*2*a), 255-Math.floor(256*(2-2*a)), 0)
            // } else {
                //     fill(255-Math.floor(256*2*a), Math.floor(256*2*a), 0)
                // }
        if (colorAct) {
            var a = Math.floor(256*this.activityValue/model_params.max_act)
            fill(a, 255-a, 0)
        } else {
            fill(color(COLOR[this.cellType]));
        }
        rect(dx*this.i, dx*this.j, dx, dx);
        pop();
    }
}