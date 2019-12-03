class Cell {

    constructor(type, id, sites, borderSites) {
        this.type = type
        this.id = id;
        // this.area = sites.length;
        this.targetArea = model_params.A[type];
        // this.perimeter = borderSites.length;
        this.sites = sites;
        this.borderSites = borderSites;
    }

    area() {
        return this.sites.length;
    }

    perimeter() {
        return this.borderSites.length;
    }

    addSite(site) {
        site.cellId = this.id;
        site.cellType = this.type;
        this.sites.push(site)
        this.borderSites.push(site)
        
        // remove sites from borderSites that are no longer on border
        var neighbors = neighborsOf(site);
        for (let v of neighbors) {
            if (v.cellId == this.id) {
                if (!v.isOnBorder()) {
                    this.removeBorderSite(site)
                }
            }
        }
    }

    removeSite(site) {
        for (let i=0; i<this.sites.length; i++) {
            if (this.sites[i].i == site.i && this.sites[i].j == site.j) {
                this.sites.splice(i, 1);
            }
        }
        this.removeBorderSite(site)

        // add sites to borderSites that are now on the border
        var neighbors = neighborsOf(site);
        for (let v of neighbors) {
            if (v.cellId == this.id) {
                if (!this.findBorderSite(v) && v.isOnBorder()) {
                    this.borderSites.push(v)
                }
            }
        }
    }

    removeBorderSite(site) {
        for (let i=0; i<this.borderSites.length; i++) {
            if (this.borderSites[i].i == site.i && this.borderSites[i].j == site.j) {
                this.borderSites.splice(i, 1);
            }
        }
    }

    findSite(site) {
        for (let i=0; i<this.sites.length; i++) {
            if (this.sites[i].i == site.i && this.sites[i].j == site.j) {
                return true;
            }
        }
        return false;
    }

    findBorderSite(site) {
        for (let i=0; i<this.borderSites.length; i++) {
            if (this.borderSites[i].i == site.i && this.borderSites[i].j == site.j) {
                return true;
            }
        }
        return false;
    }

    draw() {
        push();
        fill(color(COLOR[this.type]));
        noStroke()
        for(let site of this.sites) {
            // rect(dx*site.i, dx*site.j, dx, dx);
            site.draw()
        }
        if (drawBorders) {
            stroke(0)
            for(let site of this.borderSites) {
                if (site.i > 0) {
                    if (site.cellId != cpm.lattice_matrix[site.i-1][site.j].cellId) {
                        line(dx*site.i, dx*site.j, dx*site.i, dx*(site.j+1))
                    }
                }
                if (site.j < lattice_width-1) {
                    if (site.cellId != cpm.lattice_matrix[site.i][site.j+1].cellId) {
                        line(dx*site.i, dx*(site.j+1), dx*(site.i+1), dx*(site.j+1))
                    }
                }
                if (site.i < lattice_height-1) {
                    if (site.cellId != cpm.lattice_matrix[site.i+1][site.j].cellId) {
                        line(dx*(site.i+1), dx*site.j, dx*(site.i+1), dx*(site.j+1))
                    }
                }
                if (site.j > 0) {
                    if (site.cellId != cpm.lattice_matrix[site.i][site.j-1].cellId) {
                        line(dx*site.i, dx*site.j, dx*(site.i+1), dx*site.j)
                    }
                }
            }
        }
        pop();
    }
}