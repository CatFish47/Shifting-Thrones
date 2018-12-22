// Classes

/////////////////////////////////////////////////////////////////////////////////
// Ported from Stefan Gustavson's java implementation
// http://staffwww.itn.liu.se/~stegu/simplexnoise/simplexnoise.pdf
// Read Stefan's excellent paper for details on how this code works.
//
// Sean McCullough banksean@gmail.com

/**
 * You can pass in a random number generator object if you like.
 * It is assumed to have a random() method.
 */
var SimplexNoise = function(r) {
	if (r == undefined) r = Math;
  this.grad3 = [[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0],
                                 [1,0,1],[-1,0,1],[1,0,-1],[-1,0,-1],
                                 [0,1,1],[0,-1,1],[0,1,-1],[0,-1,-1]];
  this.p = [];
  for (var i=0; i<256; i++) {
	  this.p[i] = Math.floor(r.random()*256);
  }
  // To remove the need for index wrapping, double the permutation table length
  this.perm = [];
  for(var i=0; i<512; i++) {
		this.perm[i]=this.p[i & 255];
	}

  // A lookup table to traverse the simplex around a given point in 4D.
  // Details can be found where this table is used, in the 4D noise method.
  this.simplex = [
    [0,1,2,3],[0,1,3,2],[0,0,0,0],[0,2,3,1],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,2,3,0],
    [0,2,1,3],[0,0,0,0],[0,3,1,2],[0,3,2,1],[0,0,0,0],[0,0,0,0],[0,0,0,0],[1,3,2,0],
    [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
    [1,2,0,3],[0,0,0,0],[1,3,0,2],[0,0,0,0],[0,0,0,0],[0,0,0,0],[2,3,0,1],[2,3,1,0],
    [1,0,2,3],[1,0,3,2],[0,0,0,0],[0,0,0,0],[0,0,0,0],[2,0,3,1],[0,0,0,0],[2,1,3,0],
    [0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0],
    [2,0,1,3],[0,0,0,0],[0,0,0,0],[0,0,0,0],[3,0,1,2],[3,0,2,1],[0,0,0,0],[3,1,2,0],
    [2,1,0,3],[0,0,0,0],[0,0,0,0],[0,0,0,0],[3,1,0,2],[0,0,0,0],[3,2,0,1],[3,2,1,0]];
};

SimplexNoise.prototype.dot = function(g, x, y) {
	return g[0]*x + g[1]*y;
};

SimplexNoise.prototype.findTerrain = function(xin, yin) {
  var n0, n1, n2; // Noise contributions from the three corners
  // Skew the input space to determine which simplex cell we're in
  var F2 = 0.5*(Math.sqrt(3.0)-1.0);
  var s = (xin+yin)*F2; // Hairy factor for 2D
  var i = Math.floor(xin+s);
  var j = Math.floor(yin+s);
  var G2 = (3.0-Math.sqrt(3.0))/6.0;
  var t = (i+j)*G2;
  var X0 = i-t; // Unskew the cell origin back to (x,y) space
  var Y0 = j-t;
  var x0 = xin-X0; // The x,y distances from the cell origin
  var y0 = yin-Y0;
  // For the 2D case, the simplex shape is an equilateral triangle.
  // Determine which simplex we are in.
  var i1, j1; // Offsets for second (middle) corner of simplex in (i,j) coords
  if(x0>y0) {i1=1; j1=0;} // lower triangle, XY order: (0,0)->(1,0)->(1,1)
  else {i1=0; j1=1;}      // upper triangle, YX order: (0,0)->(0,1)->(1,1)
  // A step of (1,0) in (i,j) means a step of (1-c,-c) in (x,y), and
  // a step of (0,1) in (i,j) means a step of (-c,1-c) in (x,y), where
  // c = (3-sqrt(3))/6
  var x1 = x0 - i1 + G2; // Offsets for middle corner in (x,y) unskewed coords
  var y1 = y0 - j1 + G2;
  var x2 = x0 - 1.0 + 2.0 * G2; // Offsets for last corner in (x,y) unskewed coords
  var y2 = y0 - 1.0 + 2.0 * G2;
  // Work out the hashed gradient indices of the three simplex corners
  var ii = i & 255;
  var jj = j & 255;
  var gi0 = this.perm[ii+this.perm[jj]] % 12;
  var gi1 = this.perm[ii+i1+this.perm[jj+j1]] % 12;
  var gi2 = this.perm[ii+1+this.perm[jj+1]] % 12;
  // Calculate the contribution from the three corners
  var t0 = 0.5 - x0*x0-y0*y0;
  if(t0<0) n0 = 0.0;
  else {
    t0 *= t0;
    n0 = t0 * t0 * this.dot(this.grad3[gi0], x0, y0);  // (x,y) of grad3 used for 2D gradient
  }
  var t1 = 0.5 - x1*x1-y1*y1;
  if(t1<0) n1 = 0.0;
  else {
    t1 *= t1;
    n1 = t1 * t1 * this.dot(this.grad3[gi1], x1, y1);
  }
  var t2 = 0.5 - x2*x2-y2*y2;
  if(t2<0) n2 = 0.0;
  else {
    t2 *= t2;
    n2 = t2 * t2 * this.dot(this.grad3[gi2], x2, y2);
  }
  // Add contributions from each corner to get the final noise value.
  // The result is scaled to return values in the interval [-1,1].
  var result = 70.0 * (n0 + n1 + n2);

  // My code below
  var terrain = ["water", "plains", "forest", "mountain"]

  if (-1 <= result && result < 2 / terrain.length * 1 - 1) {
    return terrain[0];
  } else if (result < 2 / terrain.length * 2 - 1) {
    return terrain[1];
  } else if (result < 2 / terrain.length * 3 - 1) {
    return terrain[2];
  } else if (result <= 2 / terrain.length * 4 - 1) {
    return terrain[3];
  }
};

///////////////////////////////////////////////////////////////////////////////

// Player class
class Player { // TODO: Implement player depending on what class they choose
  constructor(type) {
    this.q = 0;
    this.r = -boardSize;
    this.image = new Image();
    this.image.src = "images/whiteCastle.png";
  }

  drawPlayer() {
    var x = tileSize * (this.q * 3/2) + scrollX;
    var y = tileSize * (this.q * Math.sqrt(3)/2 + this.r * Math.sqrt(3)) + scrollY;

    context.drawImage(this.image, x - 0.5 * tileSize, y - 0.5 * tileSize, tileSize, tileSize);
  }
}

// Card Class
class Card { // TODO: Add in cards depending on what the set is and determines what the ability is depending on an array of names
  constructor(set, name) {

  }
}

// Tile Class
class Tile {
  constructor(q, r, terrain) {
    this.q = q;
    this.r = r;
    this.terrain = terrain;
    this.resource = null;

    this.occupied = false;
    this.visible = true;
    this.seen = true;
  }

  generateResources() {
    var rand = Math.random();

    if (this.terrain == "forest") {
      if (rand < 0.33) {
        this.resource = new Wood();
      } else if (rand < 0.66) {
        this.resource = new Food();
      } else {
        this.resource = new Plants();
      }
    } else if (this.terrain == "mountain") {
      this.resource = new Metal();
    } else if (this.terrain == "plains") {
      if (rand < 0.5) {
        this.resource = new Food();
      } else {
        this.resource = new Plants();
      }
    }
  }

  mouseInsideTile() {
    var x = tileSize * (this.q * 3/2) + scrollX;
    var y = tileSize * (this.q * Math.sqrt(3)/2 + this.r * Math.sqrt(3)) + scrollY;

    if (mouseX < x - tileSize || mouseX > x + tileSize || mouseY > y + tileSize * Math.sqrt(3)/2 || mouseY < y - tileSize * Math.sqrt(3)/2) {
      return false; // If the mouse is not in the hexagon's bounding box, return false
    } else {
      var tempX = mouseX; // Make temporary mouse vars
      var tempY = mouseY;

      tempX = Math.abs(tempX - x) + x; // Make sure that x and y positions are in the positive quadrant
      tempY = Math.abs(tempY - y) + y;

      if (tempY - (y + tileSize * Math.sqrt(3)/2) < -1 * Math.sqrt(3) * (tempX - (x + tileSize * 1/2))) {
        return true; // Formula could use some improvement
      } else {
        return false;
      }
    }
  }

  drawTile() {
    // Translate axial coordinates to Cartesian coordinates
    var x = tileSize * (this.q * 3/2) + scrollX;
    var y = tileSize * (this.q * Math.sqrt(3)/2 + this.r * Math.sqrt(3)) + scrollY;

    // Replace with hex images
    context.beginPath();
    context.moveTo(x - tileSize * 1/2, y + tileSize * Math.sqrt(3)/2);
    context.lineTo(x + tileSize * 1/2, y + tileSize * Math.sqrt(3)/2);
    context.lineTo(x + tileSize, y);
    context.lineTo(x + tileSize * 1/2, y - tileSize * Math.sqrt(3)/2);
    context.lineTo(x - tileSize * 1/2, y - tileSize * Math.sqrt(3)/2);
    context.lineTo(x - tileSize, y);
    context.lineTo(x - tileSize * 1/2, y + tileSize * Math.sqrt(3)/2);
    context.stroke();
    context.closePath();

    if (this.terrain == "mountain") {
      context.fillStyle = "#867e70";
    } else if (this.terrain == "water") {
      context.fillStyle = "#1e5878";
    } else if (this.terrain == "plains") {
      context.fillStyle = "#b5902f";
    } else if (this.terrain == "forest") {
      context.fillStyle = "#415241";
    }

    context.fill();

    if (this.resource) {
      context.drawImage(this.resource.image, x - tileSize * 0.5, y - tileSize * 0.5, tileSize, tileSize);
    }
  }

  checkEdge() {
    if (this.q == -boardSize || this.q == boardSize) {
      return true;
    } else if (this.r == boardSize || this.r == -boardSize) {
      return true;
    } else if (this.r + this.q == boardSize || this.r + this.q == -boardSize) {
      return true;
    } else {
      return false;
    }
  }
}

// Resources Class
// TODO: Add a static function for gathering cards and for cooldowns for gathering
class Metal {
  constructor() {
    this.image = new Image();
    this.image.src = "images/metal.png";
    this.name = "metal";
  }
}
class Plants {
  constructor() {
    this.image = new Image();
    this.image.src = "images/plants.png";
    this.name = "plants";
  }
}
class Food {
  constructor() {
    this.image = new Image();
    this.image.src = "images/food.png";
    this.name = "food";
  }
}
class Wood {
  constructor() {
    this.image = new Image();
    this.image.src = "images/wood.png";
    this.name = "wood";
  }
}

// General variables
var food = [
  {}
]; // JSON for food cards -- generally for troops
var plants = []; // JSON for plant cards -- generally for healing
var metal = []; // JSON for metal cards -- generally for tools and traps
var wood = []; // JSON for wood cards -- generally for traps and ranged weapons

// Game variables
var turn = 0;
var won = false;
var board = []; // A list of all the tiles in the game
var zoom = 1.00; // A percentage
var scrollX = 0; // How far away the user's "sight" is from the center
var scrollY = 0;
var oSize = 30;
var tileSize = oSize;
var boardSize = 10;
var resourceFreq = 0.13;
var players = [];

// Viewing Variables
var scrollMode = false;
var dragMode = false;
var mouseX = 0;
var mouseY = 0;

// Canvas Variables
var $canvas = document.querySelector('canvas');
var context = $canvas.getContext("2d");
$canvas.width = 1200;
$canvas.height = 700;

document.addEventListener("wheel", function(e) {
  var zoomScale = 0.15;

  if (e.deltaY < 0) { // Scrollin Up
    if (zoom <= 3) {
      zoom += zoomScale;
      tileSize = oSize * zoom;
    } else {
      zoom = 3;
    }
  }
  if (e.deltaY > 0) { // Scrollin Down
    if (zoom >= 0.5) {
      zoom -= zoomScale;
      tileSize = oSize * zoom;
    } else {
      zoom = 0.5;
    }
  }

  document.getElementById("zoom-info").innerHTML = `x${Math.round(zoom * 10000) / 100}%`;

  // Add zoom to mouse
})

document.addEventListener("keydown", function(e) {
  if (e.keyCode == 17) { // If control key is pressed
    scrollMode = true;
  }
})

document.addEventListener("keyup", function(e) {
  if (e.keyCode == 17) {
    scrollMode = false;
  }

  if (e.keyCode == 32) {
    destroyBorders();
  }
})

document.addEventListener("mousedown", function(e) {
  if (scrollMode) {
    dragMode = true;
  }
})

document.addEventListener("mouseup", function(e) {
  dragMode = false;
})

document.addEventListener("mousemove", function(e) { // Panning
  var xMaxScale = boardSize * 100;
  var yMaxScale = boardSize * 130;

  if (dragMode) {
    if (scrollX + e.movementX > xMaxScale) {
      scrollX = xMaxScale;
    } else if (scrollX < -xMaxScale) {
      scrollX = -xMaxScale;
    } else {
      scrollX += e.movementX;
    }

    if (scrollY + e.movementY > yMaxScale) {
      scrollY = yMaxScale;
    } else if (scrollY + e.movementY < -yMaxScale) {
      scrollY = -yMaxScale;
    } else {
      scrollY += e.movementY;
    }
  }
})

document.addEventListener("mousemove", function(e) { // Tile identifying
  if (zoom > 0.9) {
    for (var i = 0; i < board.length; i++) {
      if (board[i].mouseInsideTile()) {
        document.getElementById("tile-info").innerHTML = `(${board[i].q}, ${board[i].r}): ${board[i].terrain}`;

        if (board[i].resource) {
          document.getElementById("tile-info").innerHTML += `, ${board[i].resource.name}`;
        }
        break;
      }
    }
  } else {
    document.getElementById("tile-info").innerHTML = '';
  }
})

document.getElementById("shifting-thrones").addEventListener("mousemove", function(e) { // Getting mouse positions
  mouseX = e.clientX;
  mouseY = e.clientY;
})

function drawBoard() {
  for (var i = 0; i < board.length; i++) {
    board[i].drawTile();
  }

  for (var i = 0; i < players.length; i++) {
    players[i].drawPlayer();
  }
}

function setupBoard(size) {
  var q = -size;
  var r = 0;

  while (true) {
    while (true) {
      board.push(new Tile(q, r, terrainGen.findTerrain(q, r)));

      if (q <= 0 && r == size) {
        break;
      } else if (q > 0 && r == size - q) {
        break;
      } else {
        r++;
      }
    }

    if (q == size && r == 0) {
      break;
    } else if (q < 0) {
      q++;
      r = -size - q;
    } else if (q >= 0) {
      q++;
      r = -size;
    }
  }

  for (var i = 0; i < board.length; i++) {
    var rand = Math.random();
    if (rand < resourceFreq) { // The chance for each tile to generate a resource
      board[i].generateResources();
    }
  }
}

function destroyBorders() {
  var count = 0;

  while (true) {
    if (board.length <= count) {
      break;
    }

    if (checkEdge(board[count])) {
      board.splice(count, 1);
      count -= 1;
    }

    count++;
  }

  boardSize--;
}

function init() {
  setupBoard(boardSize);
  drawBoard(board);
  window.requestAnimationFrame(animate);
}

function update() { // TODO: Add in selection of tiles
  // Animate stuff
}

function draw() {
  context.clearRect(0, 0, $canvas.width, $canvas.height);
  drawBoard();

  return 0;
}

function animate() {
  update();

  if (draw() == 0) {
    if (!won) {
      window.requestAnimationFrame(animate);
    } else {
      gameOverScreen();
    }
  }
}

var terrainGen = new SimplexNoise(Math);

players.push(new Player('white'))

init();
