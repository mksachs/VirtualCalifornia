var Popup = Class.$extend({

__init__ : function(pu_dom_id, callback_fn) {
    _this_ = this;
    
    this.pu_dom_id = pu_dom_id;
    this.pu_dom_jq = $(this.pu_dom_id);
    
    this.fade_duration = 100;
    this.popup_active = false;
    
    this.pu_dom_jq.append( "<div class=\"popup_window\"></div>" );
    this.pu_dom_jq.children(".popup_window").append( this.pu_dom_jq.children("ul") );
    this.pu_dom_jq.children("ul").remove()
    
    this.pu_dom_jq.mouseenter({_this_:_this_}, _this_.popupHovered);
    this.pu_dom_jq.mouseleave({_this_:_this_}, _this_.popupUnHovered);
    this.pu_dom_jq.click({_this_:_this_}, _this_.popupClicked);
    
    $("html").click({_this_pu_:this}, this.clickOutside);
    $(window).resize({_this_:this}, this.resizeWindow);
    
    this.setSelectedItem();
},

setSelectedItem : function(event) {
    // the Popup object
    if ( event == undefined ) {
        _this_ = this
    } else {
        _this_ = event.data._this_
    }
    
   _this_.pu_dom_jq.append(_this_.pu_dom_jq.find(".selected").text());

},

popupHovered : function(event) {
    $(event.delegateTarget).css("cursor","pointer");
},

popupUnHovered : function(event) {
    $(event.delegateTarget).css("cursor","auto");
},

popupClicked : function(event) {
    event.stopPropagation();
    _this_ = event.data._this_
    
    // clone the menu item
    subitem = $(event.delegateTarget).children(".popup_window").clone(true, true).appendTo("#uber");
    
    subitem.fadeIn(_this_.fade_duration);
    
    subitem.css({"position":"absolute"});
    
    _this_.positionPopup(event, subitem);
        
    _this_.popup_active = true;
},

positionPopup : function(event, subitem) {
    _this_ = event.data._this_
    
    offset = _this_.pu_dom_jq.offset();
    
    ////console.debug([subitem.outerWidth(), offset.left, subitem.outerWidth() + offset.left, $(window).width()]);
    
    if ( subitem.outerWidth() + offset.left < $(window).width() ) {
        subitem.offset(
            {   top:
                        offset.top +
                        _this_.pu_dom_jq.outerHeight(),
                left:
                        offset.left
            });
    } else {
        subitem.offset(
            {   top:
                        offset.top +
                        _this_.pu_dom_jq.outerHeight(),
                left:
                        offset.left - (subitem.outerWidth() - _this_.pu_dom_jq.outerWidth())
            });
    }

},

hidePopup : function(event) {
    _this_ = event.data._this_;
    $("#uber").children(".popup_window").fadeOut(_this_.fade_duration, function() {
        $("#uber").children(".popup_window").remove();
    });
    _this_.popup_active = false;
},

clickOutside : function(event) {
    _this_ = event.data._this_;
    //console.debug("popup");
    if (_this_.popup_active) {
        // hide the current menu item
        _this_.hidePopup(event);
    }
},

resizeWindow : function(event) {
    _this_ = event.data._this_;
    if (_this_.popup_active) {
        _this_.positionPopup(event, $("#uber").children(".popup_window"));
    }
}

});