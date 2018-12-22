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

SimplexNoise.prototype.noise = function(xin, yin) {
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
  return 70.0 * (n0 + n1 + n2);
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
}

// Card Class
class Card { // TODO: Add in cards depending on what the set is and determines what the ability is depending on an array of names
  constructor(set, name) {

  }
}

// Resources Class
// TODO: Add a static function for gathering cards and for cooldowns for gathering
class Metal {
  constructor() {
    this.image = new Image();
    this.image.src = "images/metal.png";
  }
}
class Plants {
  constructor() {
    this.image = new Image();
    this.image.src = "images/plants.png";
  }
}
class Food {
  constructor() {
    this.image = new Image();
    this.image.src = "images/food.png";
  }
}
class Wood {
  constructor() {
    this.image = new Image();
    this.image.src = "images/wood.png";
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
var boardSize = 15;
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

document.addEventListener("mousemove", function(e) {
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

document.getElementById("shifting-thrones").addEventListener("mousemove", function(e) {
  mouseX = e.clientX;
  mouseY = e.clientY;
})

function drawTile(q, r, size, info) { // Info is in [terrain, resources, list of traps, (Add more)]
  // Translate axial coordinates to Cartesian coordinates
  var x = size * (q * 3/2) + scrollX;
  var y = size * (q * Math.sqrt(3)/2 + r * Math.sqrt(3)) + scrollY;
  var terrain = info[0];
  var resource = info[1];

  // Replace with hex images
  context.beginPath();
  context.moveTo(x - size * 1/2, y + size * Math.sqrt(3)/2);
  context.lineTo(x + size * 1/2, y + size * Math.sqrt(3)/2);
  context.lineTo(x + size, y);
  context.lineTo(x + size * 1/2, y - size * Math.sqrt(3)/2);
  context.lineTo(x - size * 1/2, y - size * Math.sqrt(3)/2);
  context.lineTo(x - size, y);
  context.lineTo(x - size * 1/2, y + size * Math.sqrt(3)/2);
  context.stroke();
  context.closePath();

  if (terrain == "mountain") {
    context.fillStyle = "#867e70";
  } else if (terrain == "water") {
    context.fillStyle = "#1e5878";
  } else if (terrain == "plains") {
    context.fillStyle = "#b5902f";
  } else if (terrain == "forest") {
    context.fillStyle = "#415241";
  }

  context.fill();

  if (resource) {
    context.drawImage(resource.image, x - size * 0.5, y - size * 0.5, size, size);
  }
}

function drawBoard(size) {
  // For bla bla on board, drawHex

  for (var i = 0; i < board.length; i++) {
    var q = board[i][0];
    var r = board[i][1];
    var terrain = board[i][2];
    var resource = board[i][3];

    drawTile(q, r, tileSize, [terrain, resource]);
  }
}

function drawPlayers() {
  for (var i = 0; i < players.length; i++) {
    var x = tileSize * (players[i].q * 3/2) + scrollX;
    var y = tileSize * (players[i].q * Math.sqrt(3)/2 + players[i].r * Math.sqrt(3)) + scrollY;

    context.drawImage(players[i].image, x - 0.5 * tileSize, y - 0.5 * tileSize, tileSize, tileSize);
  }
}

function setupBoard(size) {
  var q = -size;
  var r = 0;

  while (true) {
    while (true) {
      board.push([q, r, findTerrain(noiseGen.noise(q, r))]);

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
      board[i].push(generateResources(board[i]));
    } else {
      board[i].push(null);
    }
  }
}

function findTerrain(noise) {
  var terrain = ["water", "plains", "forest", "mountain"]

  if (-1 <= noise && noise < 2 / terrain.length * 1 - 1) {
    return terrain[0];
  } else if (noise < 2 / terrain.length * 2 - 1) {
    return terrain[1];
  } else if (noise < 2 / terrain.length * 3 - 1) {
    return terrain[2];
  } else if (noise <= 2 / terrain.length * 4 - 1) {
    return terrain[3];
  }
}

function generateResources(tile) { // TODO: Implement a generation of resources depending on terrain
  var terrain = tile[2];
  var rand = Math.random();

  if (terrain == "water") {
    return null;
  } else if (terrain == "forest") {
    if (rand < 0.33) {
      return new Wood();
    } else if (rand < 0.66) {
      return new Food();
    } else {
      return new Plants();
    }
  } else if (terrain == "mountain") {
    return new Metal();
  } else if (terrain == "plains") {
    if (rand < 0.5) {
      return new Food();
    } else {
      return new Plants();
    }
  }
}

function destroyBorders() {
  var count = 0;

  while (true) {
    if (board.length < count) {
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

function checkEdge(tile) {
  if (tile) {
    var q = tile[0];
    var r = tile[1];

    if (q == -boardSize || q == boardSize) {
      return true;
    } else if (r == boardSize || r == -boardSize) {
      return true;
    } else if (r + q == boardSize || r + q == -boardSize) {
      return true;
    } else {
      return false;
    }
  }

  return false;
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
  drawBoard(tileSize);
  drawPlayers();

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

var noiseGen = new SimplexNoise(Math);

players.push(new Player('white'))

init();
