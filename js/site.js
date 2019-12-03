class Site {

    constructor(i, j, cellId, cellType, activityValue) {
        this.i = i;
        this.j = j;
        this.cellId = cellId;
        this.cellType = cellType; // integer
        this.activityValue = activityValue;
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
        fill(color(COLOR[this.cellType]));
        rect(dx*this.i, dx*this.j, dx, dx);
        pop();
    }
}