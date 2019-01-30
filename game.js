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
  constructor(type, id) {
    this.q = 0;
    this.r = -boardSize;
    this.image = new Image();
    this.image.src = "images/whiteCastle.png";

    this.turn = 0;

    this.maxHP = 10;
    this.hp = this.maxHP;

    this.maxMoves = 5;
    this.moves = this.maxMoves;

    this.selected = false; // When clicked on, true

    this.cards = [];
    this.traps = [];
    this.units = [];

    this.you;
    this.onResource;

    if (id == playerId) {
      this.you = true;
    } else {
      this.you = false;
    }
  }

  draw() {
    var x = boardCenter.x + tileSize * (this.q * 3/2) + scrollX;
    var y = boardCenter.y + tileSize * (this.q * Math.sqrt(3)/2 + this.r * Math.sqrt(3)) + scrollY;

    context.drawImage(this.image, x - 0.5 * tileSize, y - 0.5 * tileSize, tileSize, tileSize);
  }

  update() {
    var positions = this.calcPossibleMovement();

    for (var i = 0; i < positions.length; i++) {
      //this.selected ? positions[i].canMoveOn = true : positions[i].canMoveOn = false;
      if (this.selected) {
        positions[i].canMoveOn = true;
      }
    }

    if (findTile(this.q, this.r).resource) {
      this.onResource = true;
    } else {
      this.onResource = false;
    }
  }

  move(q, r) {
    //if (this.moves > 0) {
    this.q = q;
    this.r = r;

    for (var i = 0; i < board.length; i++) {
      board[i].canMoveOn = false;
    }

    this.moves--;
    //}
  }

  calcPossibleMovement() {
    var positions = [];
    var movePower = 2;

    for (var i = 0; i < board.length; i++) {
      if (this.inRadius(movePower, board[i]) && board[i].terrain != "water") {
        positions.push(board[i]);
      }
    }

    return positions;
  }

  inRadius(rad, tile) {
    var q = -rad;
    var r = 0;
    var tempBoard = [];

    while (true) {
      while (true) {
        tempBoard.push(new Tile(q + this.q, r + this.r, null));

        if (q <= 0 && r == rad) {
          break; // If the col# is less than 0 and r reaches the bottom, move on
        } else if (q > 0 && r == rad - q) {
          break; // If it's more than 0, and r reaches bottom, move on
        } else {
          r++; // If not, continue with increasing the row#
        }
      }

      if (q == rad && r == 0) { // If last column and last row, end
        break;
      } else if (q < 0) { // Reset to next column
        q++;
        r = -rad - q;
      } else if (q >= 0) { // Reset to next column
        q++;
        r = -rad;
      }
    }

    for (var i = 0; i < tempBoard.length; i++) {
      if (tempBoard[i].q == tile.q && tempBoard[i].r == tile.r) {
        return true;
      }
    }

    return false;
  }

  collectResource(q, r) {
    var setName = findTile(q, r).resource.name;

    if (findTile(q, r).resource.cooldown == 0) {
      for (var i = 0; i < 3; i++) {
        this.collectCard(setName);
      }
    }

    findTile(q, r).resource.maxCooldown();

    this.moves--;

    console.log(this.cards)
  }

  collectCard(setName) {
    var set = eval(setName);

    this.cards.push(set[Math.floor(Math.random() * set.length)]);
  }
}

// Card Class
class Card {
  constructor(set, name) {
    var card;

    for (var i = 0; i < set.length; i++) {
      if (set[i].name == name) {
        card = set[i];
      }
    }

    this.name = card.name;
    this.desc = card.desc;
    this.effect = eval(card.effect);
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
    this.canMoveOn = false;
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
    var x = boardCenter.x + tileSize * (this.q * 3/2) + scrollX;
    var y = boardCenter.y + tileSize * (this.q * Math.sqrt(3)/2 + this.r * Math.sqrt(3)) + scrollY;

    if (mouseX < x - tileSize || mouseX > x + tileSize ||
      mouseY > y + tileSize * Math.sqrt(3)/2 || mouseY < y - tileSize * Math.sqrt(3)/2) {
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

  draw() {
    // Translate axial coordinates to Cartesian coordinates
    var x = boardCenter.x + tileSize * (this.q * 3/2) + scrollX;
    var y = boardCenter.y + tileSize * (this.q * Math.sqrt(3)/2 + this.r * Math.sqrt(3)) + scrollY;

    if (this.canMoveOn) {
      context.strokeStyle = "#fff";
      context.lineWidth = 5;
    } else {
      context.strokeStyle = "#000";
      context.lineWidth = 1;
    }

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

    context.lineWidth = 1;


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

// Trap class
class Trap {
  constructor(q, r, range, dmg, image) {
    this.q = q;
    this.r = r;
    this.range = range;
    this.dmg = dmg;
    this.image = new Image();
    image.src = image;
  }

  inRange(player) {
    var q = -this.range;
    var r = 0;
    var tempBoard = [];

    while (true) {
      while (true) {
        tempBoard.push(new Tile(q + this.q, r + this.r, null));

        if (q <= 0 && r == this.range) {
          break; // If the col# is less than 0 and r reaches the bottom, move on
        } else if (q > 0 && r == this.range - q) {
          break; // If it's more than 0, and r reaches bottom, move on
        } else {
          r++; // If not, continue with increasing the row#
        }
      }

      if (q == this.range && r == 0) { // If last column and last row, end
        break;
      } else if (q < 0) { // Reset to next column
        q++;
        r = -this.range - q;
      } else if (q >= 0) { // Reset to next column
        q++;
        r = -this.range;
      }
    }

    for (var i = 0; i < tempBoard.length; i++) {
      if (tempBoard[i].q == player.q && tempBoard[i].r == player.r) {
        player.hp -= this.dmg;
        this.destroy();
      }
    }
  }

  destroy() {
    for (var i = 0; i < players.length; i++) {
      for (var j = 0; i < players[i].traps.length; j++) {
        if (players[i].traps[j].q == this.q && players[i].traps[j].r == this.r) {
          players[i].traps[j].pop(); // NOTE: This code might be wrong
        }
      }
    }
  }
}

// Units class
/*class Unit {
constructor(q, r, dmg, ) {

}
}*/

// Attack Class
class Attack {
  constructor(q, r, dmg) {
    this.q = q;
    this.r = r;
    this.dmg = dmg;
  }

  initiate() {
    // Add code later
  }
}

// Resources Class
// TODO: Add a static function for gathering cards and for cooldowns for gathering
class Resource {
  constructor() {
    this.cooldown = 0;
  }

  maxCooldown() {
    this.cooldown = 2;
  }

  reduceCooldown() {
    if (this.cooldown > 0) {
      this.cooldown--;
    }
  }
}

class Metal extends Resource {
  constructor() {
    super();
    this.image = new Image();
    this.image.src = "images/metal.png";
    this.name = "metal";
  }
}
class Plants extends Resource {
  constructor() {
    super();
    this.image = new Image();
    this.image.src = "images/plants.png";
    this.name = "plants";
  }
}
class Food extends Resource {
  constructor() {
    super();
    this.image = new Image();
    this.image.src = "images/food.png";
    this.name = "food";
  }
}
class Wood extends Resource {
  constructor() {
    super();
    this.image = new Image();
    this.image.src = "images/wood.png";
    this.name = "wood";
  }
}

// General variables
var playerId = "lol";

var food = [
  {
    name: "Footsoldiers",
    desc: "Places a camp of footsoldiers that will deal 1 damage to any nearby castles.",
    effect: function(q, r) {

    }}
]; // JSON for food cards -- generally for troops
var plants = [
  {
    name: "Marigold",
    desc: "Heals 1 HP",
    effect: function() {
      if (this.hp != this.maxHP) {
        this.hp += 1;
      }
      this.moves--;
    }
  },
  {
    name: "Ginseng",
    desc: "Costs 1 move. Recover 2 moves.",
    effect: function() {
      if (this.moves != this.maxMoves - 1) {
        this.moves += 2;
      } else {
        this.moves = this.maxMoves;
      }
    }
  }
]; // JSON for plant cards -- generally for healing
var metal = [
  {
    name: "Mine",
    desc: "Place a proximity mine. Enemies take 2 damage when one tile away.",
    effect: function(q, r) {
      this.traps.push(new Trap(q, r, 1, 2, '')); // Add an image later
      this.moves--;
    }
  }
]; // JSON for metal cards -- generally for tools and traps
var wood = [
  {
    name: "Trebuchet",
    desc: "Costs 2 moves. Deploy and launch a one-time use trebuchet at an enemy. Deals 3 damage.",
    effect: function(q, r) {
      var atk = new Attack(q, r, 3);
      atk.initiate();
      this.move -= 2;
    }
  }
]; // JSON for wood cards -- generally for traps and ranged weapons

// Game variables
var won = false;
var board = []; // A list of all the tiles in the game
var zoom = 1.00; // A percentage
var scrollX = 0; // How far away the user's "sight" is from the center
var scrollY = 0;
var oSize = 30;
var tileSize = oSize;
var boardSize = 5;
var resourceFreq = 0.2;
var players = [];
var infoUp = false;

// Button Condition Variables
var extract = false;

// Viewing Variables
var scrollMode = false;
var dragMode = false;
var mouseX = 0;
var mouseY = 0;
var cardMenuScroll = 0;
var cardHeight = 0.05;

// Canvas Variables
var $canvas = document.querySelector('canvas');
var context = $canvas.getContext("2d");
$canvas.width = screen.width * 0.8;
$canvas.height = screen.height * 0.7;

// Border variables
var topLeft = {
  x: $canvas.width * 0.2,
  y: $canvas.height * 0.05
}
var bottomRight = {
  x: $canvas.width,
  y: $canvas.height * 0.9
}
var boardCenter = {
  x: (topLeft.x + bottomRight.x) / 2,
  y: (topLeft.y + bottomRight.y) / 2
}

// Board Event Listeners
document.addEventListener("wheel", function(e) {// Zooming
  if (mouseX > topLeft.x && mouseX < bottomRight.x && mouseY > topLeft.y && mouseY < bottomRight.y) {
    var zoomScale = 0.1;

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
  }
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
  //if (scrollMode) {
  dragMode = true;
  //}
})

document.addEventListener("mouseup", function(e) {
  dragMode = false;
})

document.addEventListener("mousemove", function(e) { // Panning
  var xMaxScale = 10000;
  var yMaxScale = 10000;

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

/*document.addEventListener("mousemove", function(e) { // Tile identifying
if (zoom > 0.9) {
for (var i = 0; i < board.length; i++) {
if (board[i].mouseInsideTile()) {
var text = `(${board[i].q}, ${board[i].r}): ${board[i].terrain}`;

if (board[i].resource) {
text += `, ${board[i].resource.name}`;
}

context.fillStyle = "#fff"; // Info bar
context.fillRect($canvas.width * 0.84, $canvas.height * 0.53, $canvas.width * 0.15, $canvas.height * 0.35);

context.fillStyle = "#000";
context.fillText(text, $canvas.width * 0.85, $canvas.height * 0.54);

break;
}
}
}
})*/

document.addEventListener("click", function(e) { // Moving the castles
  var clicked = false;

  if (mouseX > topLeft.x && mouseX < bottomRight.x && mouseY > topLeft.y && mouseY < bottomRight.y) {
    var you = getYou();

    for (var j = 0; j < board.length; j++) {
      if (you.selected && board[j].mouseInsideTile()) {
        if (!board[j].canMoveOn) {
          clicked = true;
          break;
        }
        you.move(board[j].q, board[j].r);
        you.selected = false;
        clicked = true;
        break;
      }
    }

    if (!clicked) {
      for (var i = 0; i < board.length; i++) { // if the player clicked on a tile they want to move on
        var end = false;

        for (var j = 0; j < players.length; j++) {
        if (players[j].q == board[i].q && players[j].r == board[i].r && board[i].mouseInsideTile()) {
          players[j].selected = true;
          end = true;
          break;
        } else {
          players[j].selected = false;
        }
      }

      if (end) {break;}
    }
  }
}
})

document.getElementById("shifting-thrones").addEventListener("mousemove", function(e) {
  mouseX = e.clientX;
  mouseY = e.clientY;
})
// Mouse positions

document.addEventListener("click", function(e) { // For extract button
  if (extract &&
    mouseX > $canvas.width * 0.02 && mouseX < $canvas.width * 0.22 &&
    mouseY > $canvas.height * 0.92 && mouseY < $canvas.height * 0.98) {

      var you = getYou();
      you.collectResource(you.q, you.r);
    }
  })

document.addEventListener("wheel", function(e) {// Card Menu Scrolling
  if (mouseX > 0 && mouseX < topLeft.x && mouseY > topLeft.y && mouseY < bottomRight.y) {
    if (e.deltaY < 0) { // Scrollin Up
      if (getYou().cards.length > 17 &&
          cardMenuScroll > -(getYou().cards.length - 17) * $canvas.height * 0.05) {
        cardMenuScroll -= $canvas.height * 0.05 / 8;
      }
    } else {
      if (cardMenuScroll < 0) {
        cardMenuScroll += $canvas.height * 0.05 / 8;
      }
    }
  }
})

function getYou() {
  for (var i in players) {
    if (players[i].you) {
      return players[i];
    }
  }
}

function drawInfo(li) {

}

function drawBoard() {
  for (var i = 0; i < board.length; i++) {
    board[i].draw();
  }

  for (var i = 0; i < players.length; i++) {
    players[i].draw();
  }
}

function drawCards() { // Draw the card meny on the side (left)
  var cards = getYou().cards;

  context.font = "15px Georgia";

  for (var i in cards) {
    context.fillStyle = "#d2b38c";
    context.fillRect(0, cardMenuScroll + $canvas.height * 0.05 + i * $canvas.height * cardHeight,
      $canvas.width * 0.2, $canvas.height * (cardHeight - 0.001));
    context.fillStyle = "#000";
    context.fillText(`${cards[i].name}`, $canvas.width * 0.005,
      cardMenuScroll + $canvas.height * 0.09 + i * $canvas.height * cardHeight);
  }
}

function drawActions() {
  // Buttons (in order):
  // Extract Resource | Play Selected Card | *TBD* | End Turn
  var resourceCooldown;

  var you = getYou(); // Is the player on the resource?

  if (you.onResource) {
    extract = true;
    resourceCooldown = findTile(you.q, you.r).resource.cooldown;
  } else {
    extract = false;
  }

  if (extract) {
    context.fillStyle = "#fff";
    context.fillRect($canvas.width * 0.02,
      $canvas.height * 0.92, $canvas.width * 0.2, $canvas.height * 0.06);
    context.fillStyle = "#000";
    context.font = "20px Georgia";

    if (resourceCooldown == 0) {
      context.fillText("Extract Resource", $canvas.width * 0.02, $canvas.height * 0.97);
    } else {
      context.fillText(`Cooldown: ${resourceCooldown}`, $canvas.width * 0.02, $canvas.height * 0.97);
    }
  }
}

function setupBoard(size) {
  var q = -size; // Represents columns
  var r = 0; // Represents diagonal rows

  while (true) {
    while (true) {
      board.push(new Tile(q, r, terrainGen.findTerrain(q + 100, r + 100))); // Add a new tile

      if (q <= 0 && r == size) {
        break; // If the col# is less than 0 and r reaches the bottom, move on
      } else if (q > 0 && r == size - q) {
        break; // If it's more than 0, and r reaches bottom, move on
      } else {
        r++; // If not, continue with increasing the row#
      }
    }

    if (q == size && r == 0) { // If last column and last row, end
      break;
    } else if (q < 0) { // Reset to next column
      q++;
      r = -size - q;
    } else if (q >= 0) { // Reset to next column
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

    if (board[count].checkEdge()) {
      board.splice(count, 1);
      count -= 1;
    }

    count++;
  }

  boardSize--;
}

function findTile(q, r) {
  for (var i = 0; i < board.length; i++) {
    if (board[i].q == q && board[i].r == r) {
      return board[i];
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
  for (var i = 0; i < players.length; i++) {
    players[i].update();
  }
}

function draw() {
  context.clearRect(0, 0, $canvas.width, $canvas.height);

  context.fillStyle = "#99f"; // board
  context.fillRect(topLeft.x, topLeft.y, $canvas.width - topLeft.x,
    $canvas.height * (bottomRight.y - topLeft.y));
  // Replace fillRect with an image of... a void? IDK
  drawBoard();

  context.fillStyle = "#000"; // Card menu
  context.fillRect(0, $canvas.height * 0.05, $canvas.width * 0.2, $canvas.height * 0.85);
  // Replace fillRect with suitable image
  drawCards();

  context.fillStyle = "#9f9"; // Game bar
  context.fillRect(0, 0, $canvas.width, $canvas.height * 0.05);
  // Replace fillRect with suitable image

  context.fillStyle = "#f99"; // Action bar
  context.fillRect(0, $canvas.height * 0.9, $canvas.width, $canvas.height * 0.1);
  // Replace fillRect with suitable image
  drawActions();

  return 0;
}

function animate() {
  // Add a state machine in update for what part of the game the player is looking at
  // i.e. The selection screen, home page, game board, etc.
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

//var testTile = new Tile(0, 0, "plains");

players.push(new Player('white', 'lol'));

players[0].collectCard("metal");

init();
