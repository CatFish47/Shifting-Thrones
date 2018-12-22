var mongoose = require('mongoose');

var model = mongoose.model('player', new mongoose.Schema({
    username: String,
    id: String,
    data: {type: mongoose.Schema.Types.Mixed},
    index: Number,
    state: {type: String, enum: ['template', 'blank'], default: 'initial'},
    category: {type: String},
    project: {type: mongoose.Schema.ObjectId, ref: 'Projects', required: true}
}));


exports.getModel = function() {
	return model;
}
