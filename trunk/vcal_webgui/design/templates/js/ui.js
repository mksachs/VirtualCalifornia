Object.size = function(obj) {
    var size = 0, key;
    for (key in obj) {
        if (obj.hasOwnProperty(key)) size++;
    }
    return size;
};

var UIEle = Class.$extend({
    __classvars__ : {
        UI_ELEMENT_TABLE : [],
        UI_POPUP_TABLE : [],
        UI_OPEN_POPUP : undefined
    },
    
    __init__ : function(dom) {
        this.dom_jq = $(dom);
        UIEle.UI_ELEMENT_TABLE.push(this);
    },

    UIEleHovered : function(event) {
        //$(event.delegateTarget).css("cursor","pointer");
        ////console.debug(event.data.$this.testvar)
    },
    
    UIEleUnhovered : function(event) {
        //$(event.delegateTarget).css("cursor","auto");
        ////console.debug(event.data.$this.testvar)
    },
    
    UIEleClicked : function(event) {
        /*$this = event.data.$this;
        if ( $this instanceof UIMenubarEle ) {
            //console.debug("menu bar");
        }
        else if ( $this instanceof UIPopupEle ) {
            //console.debug("popup");
        }*/
        ////console.debug(event.data.$this.testvar)
    }

});

var UIMenubar = UIEle.$extend({
    __init__ : function(dom, data_target) {
        this.$super(dom);
        
        this.menus = [];
        
        mb = this;
        this.dom_jq.children(".menu_bar_item").each(function( index ) {
            mb.menus.push(new UIMenubarPopup(this, mb));
        });
        
        this.data_target = data_target;
    }
});

var UIPopup = UIEle.$extend({
    __init__ : function(dom) {
        this.$super(dom);
        
        this.dom_jq.append( "<div class=\"popup_window\"></div>" );
        this.dom_jq.children(".popup_window").append( this.dom_jq.children("ul") );
        this.dom_jq.children("ul").remove();
        
        this.dom_jq.mouseenter({$this:this}, this.UIEleHovered);
        this.dom_jq.mouseleave({$this:this}, this.UIEleUnhovered);
        this.dom_jq.click({$this:this}, this.UIEleClicked);
        $("html").click(this.clickOutside);
        $(window).resize(this.resizeWindow);
        
        // initilaize all of the sub menu items
        pu = this;
        this.dom_jq.find("li").each(function(index){
            if (!$(this).hasClass("disabled") && !$(this).hasClass("active") ) {
                pu.setSubStatus($(this), "on");
            }
        });
        
        this.fade_duration = 50;
        
        UIEle.UI_POPUP_TABLE.push(this);
    },
    
    setSubStatus : function(sub_item, status) {
        switch ( status ) {
            case "on":
                sub_item.mouseenter({$this:this}, this.UIEleSubHovered);
                sub_item.mouseleave({$this:this}, this.UIEleSubUnhovered);
                sub_item.click({$this:this}, this.UIEleSubClicked);
                if ( sub_item.hasClass("disabled") )
                    sub_item.removeClass("disabled");
                if ( sub_item.hasClass("active") )
                    sub_item.removeClass("active");
                break;
            case "disabled":
                if ( sub_item.hasClass("active") )
                    sub_item.removeClass("active");
                if ( !sub_item.hasClass("disabled") )
                    sub_item.addClass("disabled");
                sub_item.off("click");
                sub_item.off("mouseenter");
                sub_item.off("mouseleave");
                break;
            case "active":
                if ( sub_item.hasClass("disabled") )
                    sub_item.removeClass("disabled");
                if ( sub_item.hasClass("active") )
                    sub_item.addClass("active");
                sub_item.off("click");
                sub_item.off("mouseenter");
                sub_item.off("mouseleave");
                break;
        }
        
    },
    
    UIEleSubHovered : function(event) {
        $(event.delegateTarget).css("cursor","pointer");
        $(event.delegateTarget).addClass("hovered");
    },

    UIEleSubUnhovered : function(event) {
        $(event.delegateTarget).css("cursor","auto");
        $(event.delegateTarget).removeClass("hovered");
    },
    
    UIEleSubClicked : function(event) {
        event.stopPropagation();
        
        $this = event.data.$this;
        
        $(event.delegateTarget).removeClass("hovered");
        setTimeout(function() {
            $(event.delegateTarget).addClass("active");
            setTimeout(function() {
                $(event.delegateTarget).removeClass("active");
                $this.hidePopup();
                $this.doAction(event);
        },100);
        },100);
        
    },
    
    clickOutside : function(event) {
        if ( UIEle.UI_OPEN_POPUP )
            UIEle.UI_OPEN_POPUP.hidePopup();
    },
    
    resizeWindow : function(event) {
        if ( UIEle.UI_OPEN_POPUP ) {
            UIEle.UI_OPEN_POPUP.positionPopup();
        }
    },
    
    UIEleHovered : function(event) {
        $this = event.data.$this;
        
        $(event.delegateTarget).css("cursor","pointer");
        if ( UIEle.UI_OPEN_POPUP != $this )
            $(event.delegateTarget).addClass("hovered");
    },
    
    UIEleUnhovered : function(event) {
        $this = event.data.$this;
        
        $(event.delegateTarget).css("cursor","auto");
        if ( UIEle.UI_OPEN_POPUP != $this )
            $(event.delegateTarget).removeClass("hovered");
    },
    
    UIEleClicked : function(event) {
        event.stopPropagation();
        
        $this = event.data.$this;
        
        if ( UIEle.UI_OPEN_POPUP ) {
            if ( UIEle.UI_OPEN_POPUP == $this ) {
                UIEle.UI_OPEN_POPUP.hidePopup();
            } else {
                UIEle.UI_OPEN_POPUP.hidePopup($this);
            }
        } else {
            $this.showPopup();
        }
    },
    
    hidePopup : function(calling_obj) {
        //$this = this;
        
        UIEle.UI_OPEN_POPUP.dom_jq.removeClass("active");
        
        UIEle.UI_OPEN_POPUP.active_popup_jq.fadeOut( UIEle.UI_OPEN_POPUP.fade_duration, function() {
            UIEle.UI_OPEN_POPUP.active_popup_jq.remove();
            UIEle.UI_OPEN_POPUP.active_popup_jq = undefined;
            UIEle.UI_OPEN_POPUP = undefined;
            if ( calling_obj )
                 calling_obj.showPopup();
        });
    },
    
    showPopup : function() {
        UIEle.UI_OPEN_POPUP = this;
        
        this.dom_jq.removeClass("hovered");
        this.dom_jq.addClass("active");
        
        // clone the menu item
        this.active_popup_jq = this.dom_jq.children(".popup_window").clone(true, true).appendTo("#uber");
        
        this.active_popup_jq.fadeIn($this.fade_duration);
        this.active_popup_jq.css({"position":"absolute"});
        
        this.positionPopup();
    },
    
    positionPopup : function() {
        if ( this.dom_jq.offset().left + this.active_popup_jq.outerWidth() >= $(window).width() ) {
            //the popup will appear off the right edge of the window
            //align the right edge of the popup with the right edge of the parent element
            left = this.dom_jq.offset().left + this.dom_jq.outerWidth() - this.active_popup_jq.outerWidth();
        } else {
            //position the popup normally
            ////align the left edge of the popup with the left edge of the parent element
            left = this.dom_jq.offset().left;
        }
        
        if ( this.dom_jq.offset().top + this.dom_jq.outerHeight() + this.active_popup_jq.outerHeight() >= $(window).height() ) {
            //the popup will appear off the bottom of the screen
            //align the bottom of the popup with the bottom of the parent element
            top = this.dom_jq.offset().top + this.dom_jq.outerHeight() - this.active_popup_jq.outerHeight();
        } else {
            //position the popup normally
            //align the top of the popup with the bottom of the parent element
            top = this.dom_jq.offset().top + this.dom_jq.outerHeight();
        }
        
        this.active_popup_jq.offset({ top: top, left: left });
    }

});

var UIMenubarPopup = UIPopup.$extend({
    __init__ : function(dom, parent_menubar) {
        this.$super(dom);
        
        this.parent_menubar = parent_menubar;
        this.name = this.dom_jq.children("p").text();
        
        //This is for testing
        /*
        action = this.dom_jq.find(".menu_bar_action");
        if (action.hasClass('dialog'))
            new UIModalDialog(action);
        */
    },
    
    doAction : function(event) {
        action = $(event.delegateTarget).children(".menu_bar_action");
        ////console.debug(action.hasClass('dialog'));
        //link = $(event.delegateTarget).children(".menu_bar_link");
        if (action.hasClass('dialog'))
            new UIModalDialog(action);
        if (action.hasClass('link'))
            document.location.href = action.children("a").attr("href");
            
            
    },
    
    showPopup : function() {
        this.$super();
        
        ////console.debug(this.name);
        if ( this.name == "Action" ) {
            if ( Object.size(this.parent_menubar.data_target.server_data.selected) > 0 ) {
                this.setSubStatus(this.active_popup_jq.find(".disabled"), "on");
            }
        }
        ////console.debug(Object.size(this.parent_menubar.data_target.server_data.selected));
    }
});

var UIButton = UIEle.$extend({
    __init__ : function(dom, callBack, callBackParams) {
        this.$super(dom);
        
        this.callBack = callBack;
        
        ////console.debug(data);
        
        this.callBackParams = callBackParams;
        
        this.dom_jq.mouseenter({$this:this}, this.UIEleHovered);
        this.dom_jq.mouseleave({$this:this}, this.UIEleUnhovered);
        this.dom_jq.click({$this:this}, this.UIEleClicked);
        
    },
    
    UIEleHovered : function(event) {
        $(event.delegateTarget).css("cursor","pointer");
        $(event.delegateTarget).addClass("hovered");
    },
    
    UIEleUnhovered : function(event) {
        $(event.delegateTarget).css("cursor","auto");
        $(event.delegateTarget).removeClass("hovered");
    },
    
    UIEleClicked : function(event) {
        $this = event.data.$this;
        
        $this.callBack($this.callBackParams);
    }
});

var UICheckbox = UIEle.$extend({
    __init__ : function(dom, callBack, callBackParams) {
        //build the checkbox. this will be sent to the UIEle init
        name = $(dom).attr("name");
        label_text = $(dom).siblings("label").text();
        
        $(dom).after("<div class=\"checkbox\"><div class=\"checkbox_box\"></div><div class=\"checkbox_label\">{{label text}}</div></div>".replace("{{label text}}",label_text));
        new_dom = $(dom).siblings(".checkbox");
        $(dom).siblings("label").remove();
        $(dom).remove();
        
        //now call the super with the new dom element.
        this.$super(new_dom);
        
        this.callBack = callBack;
        this.name = name;
        this.selected = false;
        
        this.callBackParams = callBackParams;
        
        this.dom_jq.mouseenter({$this:this}, this.UIEleHovered);
        this.dom_jq.mouseleave({$this:this}, this.UIEleUnhovered);
        this.dom_jq.click({$this:this}, this.UIEleClicked);
        
    },
    
    UIEleHovered : function(event) {
        $(event.delegateTarget).css("cursor","pointer");
        $(event.delegateTarget).addClass("hovered");
    },
    
    UIEleUnhovered : function(event) {
        $(event.delegateTarget).css("cursor","auto");
        $(event.delegateTarget).removeClass("hovered");
    },
    
    UIEleClicked : function(event) {
        $this = event.data.$this;
        
        if ( $this.selected ) {
            $(event.delegateTarget).removeClass("selected");
            $this.selected = false;
        } else {
            $(event.delegateTarget).addClass("selected");
            $this.selected = true;
        }
        
        //$this.callBack($this.callBackParams);
    }
});


var UIModalDialog = UIEle.$extend({
    __init__ : function(dom) {
        this.$super(dom);
        
        // we always have these two
        new UIButton(this.dom_jq.find( "[name='cancel']" ), this.doCancel, {$this:this});
        new UIButton(this.dom_jq.find( "[name='go']" ), this.doGo, {$this:this});
        
        //look for checkboxes
        this.dom_jq.find( "[type='checkbox']" ).each( function (index) {
            new UICheckbox(this);
        });
        
        
        $(window).resize({$this:this}, this.resizeWindow);
        
        this.fade_duration = 100;
        
        this.showDialog();
        
        },
    
    doCancel : function(params) {
        $this = params.$this;
        
        $this.hideDialog();
    },
    
    doGo : function(params) {
        $this = params.$this;
        
        //$this.hideDialog();
    },
    
    showDialog : function() {
        // create an overlay to capture all events outside the dialog box
        this.overlay = $("<div class=\"overlay\"></div>").appendTo("#uber");
        
        this.overlay.css({"position":"absolute", "display":"none"});
        
        this.overlay.offset(
            {   top:    0,
                left:   0
        });
        
        this.overlay.width($(window).width());
        this.overlay.height($("html").height());
        
        this.overlay.fadeIn($this.fade_duration);
        
        // clone the menu item
        this.active_dialog = this.dom_jq.clone(true, true).appendTo("#uber");
        
        this.active_dialog.fadeIn(this.fade_duration);
        this.active_dialog.css({"position":"absolute"});
        
        this.positionDialog();

    },
    
    hideDialog : function() {
        this.active_dialog.fadeOut(this.fade_duration);
        $this = this;
        this.overlay.fadeOut(this.fade_duration, function() {
            $this.active_dialog.remove();
            $this.overlay.remove();
        });
    },
    
    positionDialog : function() {
        this.active_dialog.offset(
            {   top:    $(window).height()/2.0 - this.active_dialog.outerHeight()/2.0,
                left:   $(window).width()/2.0 - this.active_dialog.outerWidth()/2.0
        });
    },
    
    resizeWindow : function(event) {
        $this = event.data.$this;
        
        $this.positionDialog();
    }

});

var UISingleSelectPopup = UIPopup.$extend({
    __init__ : function(dom, callBack, callBackParams) {
        this.$super(dom);
        
        this.callBackParams = callBackParams;
        
        this.setSelection(this.dom_jq.find(".selected").text());
        
        this.callBack = callBack;
        
    },
    
    setSelection : function(value) {
        selection = this.dom_jq.children("p");
        if ( selection.length == 0 ) {
            this.dom_jq.append("<p>"+value+"</p>");
        } else {
            this.dom_jq.children("p").text(value);
        }
        
        //this.dom_jq.append("<p>"+value+"</p>");
        this.callBackParams.selection = value;
    },
    
    doAction : function(event) {
        this.setSelection($(event.delegateTarget).text());
        
        this.callBack(this.callBackParams);
        
    }

});