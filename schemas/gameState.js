var mongoose = require('mongoose');

var model = mongoose.model('gameState', new mongoose.Schema({
    id: Number,
    data: {type: mongoose.Schema.Types.Mixed},
    playerCount: Number,
    players: {player1: Player, player2: Player, player3: Player, player4: Player, player5: Player, player6: Player}
}));


exports.getModel = function() {
	return model;
}
