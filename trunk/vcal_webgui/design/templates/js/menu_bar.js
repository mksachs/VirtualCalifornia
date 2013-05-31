var MenuBarEle = UIEle.$extend({

__init__ : function(dom_id) {
    this.$super(dom_id);
    
    _this_ = this;
    
    this.mb_dom_id = mb_dom_id;
    this.mb_dom_jq = $(this.mb_dom_id);
    
    this.fade_duration = 100;
    this.active_menu = null;
    this.menu_bars = new Array();
    
    this.active_menu_bar_subitem = null;
    this.active_menu_bar = null;
    
    this.in_active_subitem = false;
    
    this.active_dialog = null;
    
    this.selected_actions = new Object();
    
    this.mb_dom_jq.find(".menu_bar_item").each(function( index ) {
        // event callbacks
        $(this).mouseenter({_this_:_this_}, _this_.menuBarItemInactiveHovered);
        $(this).mouseleave({_this_:_this_}, _this_.menuBarItemInactiveUnHovered);
        $(this).click({_this_:_this_}, _this_.menuBarItemClicked);
        
        // add the divider if ness
        $(this).find("ul").each(function( index ) {
            if ( index > 0 )
                $(this).addClass("divider");
        });
        
        // store the index of the menu and put the jq obj in the menu array
        $(this).data("index", index);
        _this_.menu_bars[index] = {main:$(this), sub:null};
        
        // initilaize all of the sub menu items
        $(this).find("li").each(function(index){
            $(this).mouseenter({_this_:_this_}, _this_.menuBarItemInactiveHovered);
            $(this).mouseleave({_this_:_this_}, _this_.menuBarItemInactiveUnHovered);
            $(this).click({_this_:_this_}, _this_.menuBarSubitemClicked);
        });
        
        // initialize all of the dialogs
        $(this).find(".menu_bar_dialog").each(function(index){
            //cancel buttons always do the same thing
            cancel_button = $(this).find( "[name='cancel']" );
            cancel_button.click({_this_:_this_}, _this_.cancelButtonClicked);
        });
        
    });
    
    $("html").click({_this_mb_:this}, _this_.clickOutside);
    $(window).resize({_this_:this}, this.resizeWindow);
    
    // this is just for testing
    //this.showDialog({data:{_this_:this}, delegateTarget: $("#action_menu_select_actions")});
},

cancelButtonClicked : function(event) {
    _this_ = event.data._this_;
    
    _this_.active_dialog.fadeOut(_this_.fade_duration);
    _this_.overlay.fadeOut(_this_.fade_duration, function() {
        _this_.active_dialog.remove();
        _this_.overlay.remove();
    });

},

//Inactive hover functions
menuBarItemInactiveHovered : function(event) {
    $(event.delegateTarget).css("cursor","pointer");
    $(event.delegateTarget).addClass("hovered");
},

menuBarItemInactiveUnHovered : function(event) {
    $(event.delegateTarget).css("cursor","auto");
    $(event.delegateTarget).removeClass("hovered");
},

//Click functions
menuBarSubitemClicked : function(event) {
    event.stopPropagation();
    _this_ = event.data._this_;
    
    _this_.hideMenu(_this_.active_menu);
    
    if( $(event.delegateTarget).find(".menu_bar_dialog").length != 0 ) {
        _this_.showDialog(event);
    }

},

resizeWindow : function(event) {
    _this_ = event.data._this_;
    
    if (_this_.active_dialog != null) {
        _this_.active_dialog.offset(
            {   top:    $(window).height()/2.0 - _this_.active_dialog.outerHeight()/2.0,
                left:   $(window).width()/2.0 - _this_.active_dialog.outerWidth()/2.0
        });
        _this_.overlay.width($(window).width());
        _this_.overlay.height($(window).height());
    }

},

showDialog : function(event) {
    _this_ = event.data._this_;
    
    // create an overlay to capture all events outside the dialog box
    overlay = $("<div class=\"overlay\"></div>").appendTo("#uber");
    
    overlay.css({"position":"absolute", "display":"none"});
    
    overlay.offset(
        {   top:    0,
            left:   0
    });
    
    overlay.width($(window).width());
    overlay.height($("html").height());
    
    overlay.fadeIn(_this_.fade_duration);
    
    // clone the menu item
    dialog = $(event.delegateTarget).find(".menu_bar_dialog").clone(true, true).appendTo("#uber");
    
    dialog.fadeIn(_this_.fade_duration);
    dialog.css({"position":"absolute"});
    
    ////console.debug([$(window).width(), $(window).height()])
    
    dialog.offset(
        {   top:    $(window).height()/2.0 - dialog.outerHeight()/2.0,
            left:   $(window).width()/2.0 - dialog.outerWidth()/2.0
    });
    
    _this_.active_dialog = dialog;
    _this_.overlay = overlay;
},

menuBarItemClicked : function(event) {
    event.stopPropagation();
    _this_ = event.data._this_;
    
    // if there is an active menu hide it first. hide menu will call show menu.
    // if not just call show menu
    if (_this_.active_menu != null) {
        _this_.hideMenu(_this_.active_menu, event);
    } else {
        _this_.showMenu(event);    
    }
},

clickOutside : function(event) {
    _this_ = event.data._this_;
    //console.debug("menubar");
    if (_this_.active_menu != null) {
        // hide the current menu item
        _this_.hideMenu(_this_.active_menu, event);
    }
},

showMenu : function(event) {
    _this_ = event.data._this_;
    
    // remove the inactive hover event handlers and add the active ones for all menu items
    _this_.mb_dom_jq.find(".menu_bar_item").each(function( index ) {
        $(this).off("mouseenter");
        $(this).off("mouseleave");
        $(this).mouseenter({_this_:_this_}, _this_.menuBarItemHovered);
        $(this).mouseleave({_this_:_this_}, _this_.menuBarItemUnHovered);
    });
    
    // remove the hovered class and add the active one
    $(event.delegateTarget).removeClass("hovered");
    $(event.delegateTarget).addClass("active");
    
    // clone the menu item
    subitem = $(event.delegateTarget).children(".menu_bar_subitem").clone(true, true).appendTo("#uber");
    
    if (subitem.length == 0) {
        // no sub nav
        subitem.remove()
    } else {
        subitem.fadeIn(_this_.fade_duration);
        offset = $(event.delegateTarget).offset();
        subitem.css({"position":"absolute"});
        //subitem.width($(event.delegateTarget).children(".menu_bar_subitem").width());
        subitem.offset(
            {   top:
                        offset.top +
                        $(event.delegateTarget).height() +
                        parseFloat($(event.delegateTarget).css("padding-top").split("px")[0]) +
                        parseFloat($(event.delegateTarget).css("padding-bottom").split("px")[0]),
                left:
                        offset.left
            });
        _this_.active_menu = $(event.delegateTarget).data("index");
        this.menu_bars[_this_.active_menu].sub = subitem;
    }

},

hideMenu : function(index, event) {
    if ( event == undefined ) {
        _this_ = this;
    } else {
        _this_ = event.data._this_;
    }
    _this_.menu_bars[index].sub.fadeOut(_this_.fade_duration, function() {
        $("#uber").children(".menu_bar_subitem").remove();
        //console.debug(_this_);
        _this_.menu_bars[index].main.removeClass("active");
        _this_.menu_bars[index].sub.remove();
        _this_.menu_bars[index].sub = null;
        _this_.active_menu = null;
    
    
        // remove the active hover event handlers and add the inactive ones for all menu items
        _this_.mb_dom_jq.find(".menu_bar_item").each(function( index ) {
            $(this).off("mouseenter");
            $(this).off("mouseleave");
            $(this).mouseenter({_this_:_this_}, _this_.menuBarItemInactiveHovered);
            $(this).mouseleave({_this_:_this_}, _this_.menuBarItemInactiveUnHovered);
        });
        
        if ( event != undefined && ($(event.delegateTarget).data("index") != index)) {
            _this_.showMenu(event);
        }

    });
},

//Active hover functions
menuBarItemHovered : function(event) {
    $(event.delegateTarget).css("cursor","pointer");
    //event.stopPropagation();
    
    //_this_ = event.data._this_
    
    //$(event.delegateTarget).css("cursor","pointer");
    //$(event.delegateTarget).addClass("hovered");
    
    
    // if the menu is active and the mouse enters from the bottom do nothing
    //if ( _this_.menu_active ) {
    //    //console.debug([event.pageY , $(event.delegateTarget).offset().top+$(event.delegateTarget).height()]);
    //     if ( event.pageY < $(event.delegateTarget).offset().top+$(event.delegateTarget).height() ) {
    //        _this_.showMenu(event);
    //    }
   // }

    
},

menuBarItemUnHovered : function(event) {
     $(event.delegateTarget).css("cursor","auto");
    //event.stopPropagation();
    
    //_this_ = event.data._this_
    
    //$(event.delegateTarget).css("cursor","auto");
    
    // if the menu is not active or if it is active and we mouse off to the top left or right
    //if (!_this_.menu_active || (_this_.menu_active && event.pageY < $(event.delegateTarget).offset().top+$(event.delegateTarget).height() ) ) {
    //    $(event.delegateTarget).removeClass("hovered");
    //}
    
},


menuBarSubitemHovered : function(event) {
    _this_ = event.data._this_;
    
    _this_.in_active_subitem = true;
},

menuBarSubitemUnHovered : function(event) {
    _this_ = event.data._this_;
    
    _this_.in_active_subitem = false;
}


});