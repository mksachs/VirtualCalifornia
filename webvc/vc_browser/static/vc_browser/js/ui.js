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
    
    __init__ : function(dom, enabled) {
        this.dom_jq = $(dom);
        UIEle.UI_ELEMENT_TABLE.push(this);
        
        if ( enabled == undefined ) {
            this.enabled = true;
        } else {
            this.enabled = enabled;
        }
        
        if ( !this.enabled )
            this.dom_jq.addClass("disabled");
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
    },
    
    UIEleToggleEnabled : function(context) {
        if ( context == undefined ) {
            $this = this;
        } else {
            $this = context;
        }
        if ( $this.enabled ) {
            $this.enabled = false;
            $this.dom_jq.addClass("disabled");
        } else {
            $this.enabled = true;
            $this.dom_jq.removeClass("disabled");
        }
    },
    
    UIEleEnable : function(context) {
        if ( context == undefined ) {
            $this = this;
        } else {
            $this = context;
        }

        $this.enabled = true;
        $this.dom_jq.removeClass("disabled");
    },
    
    UIEleDisable : function(context) {
        if ( context == undefined ) {
            $this = this;
        } else {
            $this = context;
        }

        $this.enabled = false;
        $this.dom_jq.addClass("disabled");
    }


});

var UIMenubar = UIEle.$extend({
    __init__ : function(dom, data_target) {
        this.$super(dom);
        
        this.menus = [];
        
        this.data_target = data_target;
        
        mb = this;
        this.dom_jq.children(".menu_bar_item").each(function( index ) {
            mb.menus.push(new UIMenubarPopup(this, mb));
        });
    }
});

var UIPopup = UIEle.$extend({
    __init__ : function(dom, enabled) {
        this.$super(dom, enabled);
        
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
            p_left = this.dom_jq.offset().left + this.dom_jq.outerWidth() - this.active_popup_jq.outerWidth();
        } else {
            //position the popup normally
            ////align the left edge of the popup with the left edge of the parent element
            p_left = this.dom_jq.offset().left;
        }
        
        if ( this.dom_jq.offset().top + this.dom_jq.outerHeight() + this.active_popup_jq.outerHeight() >= $(window).height() ) {
            //the popup will appear off the bottom of the screen
            //align the bottom of the popup with the bottom of the parent element
            p_top = this.dom_jq.offset().top + this.dom_jq.outerHeight() - this.active_popup_jq.outerHeight();
        } else {
            //position the popup normally
            //align the top of the popup with the bottom of the parent element
            p_top = this.dom_jq.offset().top + this.dom_jq.outerHeight();
        }
        
        
        
        this.active_popup_jq.css({ "top": p_top, "left": p_left });
    }

});

var UIMenubarPopup = UIPopup.$extend({
    __init__ : function(dom, parent_menubar, callBack, callBackParams) {
        this.$super(dom);
        
        this.parent_menubar = parent_menubar;
        this.name = this.dom_jq.children("p").text();
        
        this.callBackParams = callBackParams;
        this.callBack = callBack;
        
        //This is for testing
        //action = this.dom_jq.find(".menu_bar_action");
        //if (action.hasClass('dialog'))
        //    new UIModalDialog(action, this.parent_menubar.data_target.server_data);
        
    },
    
    doAction : function(event) {
        $this = event.data.$this;
        action = $(event.delegateTarget).children(".menu_bar_action");
        ////console.debug(action.hasClass('dialog'));
        //link = $(event.delegateTarget).children(".menu_bar_link");
        if (action.hasClass('dialog'))
            new UIModalDialog(action, $this.parent_menubar.data_target.server_data);
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
    __init__ : function(dom, callBack, callBackParams, enabled) {
        this.$super(dom, enabled);
        
                
        //this.UIEleToggleActive();
        
        this.callBack = callBack;
        
        this.name = this.dom_jq.attr("name");
        
        this.callBackParams = callBackParams;
        
        this.dom_jq.mouseenter({$this:this}, this.UIEleHovered);
        this.dom_jq.mouseleave({$this:this}, this.UIEleUnhovered);
        this.dom_jq.click({$this:this}, this.UIEleClicked);
        
    },
    
    UIEleHovered : function(event) {
        $this = event.data.$this;
        
        if ( $this.enabled ) {
            $(event.delegateTarget).css("cursor","pointer");
            $(event.delegateTarget).addClass("hovered");
        }
    },
    
    UIEleUnhovered : function(event) {
        $this = event.data.$this;
        
        if ( $this.enabled ) {
            $(event.delegateTarget).css("cursor","auto");
            $(event.delegateTarget).removeClass("hovered");
        }
    },
    
    UIEleClicked : function(event) {
        $this = event.data.$this;
        
        if ( $this.enabled ) {
            $this.callBack($this.callBackParams);
        }
    }
});

var UICheckbox = UIEle.$extend({
    __init__ : function(dom, callBack, callBackParams, enabled) {
        //build the checkbox. this will be sent to the UIEle init
        name = $(dom).attr("name");
        label_text = $(dom).siblings("label").text();
        
        $(dom).after("<div class=\"checkbox\"><div class=\"checkbox_box\"></div><div class=\"checkbox_label\">{{label text}}</div></div>".replace("{{label text}}",label_text));
        new_dom = $(dom).siblings(".checkbox");
        $(dom).siblings("label").remove();
        $(dom).remove();
        
        //now call the super with the new dom element.
        this.$super(new_dom, enabled);
        
        this.callBack = callBack;
        this.name = name;
        this.selected = false;
        this.value = this.selected
        
        this.callBackParams = callBackParams;
        
        this.dom_jq.mouseenter({$this:this}, this.UIEleHovered);
        this.dom_jq.mouseleave({$this:this}, this.UIEleUnhovered);
        this.dom_jq.click({$this:this}, this.UIEleClicked);
        
        ////console.debug(name);
        
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
            $this.value = $this.selected
        } else {
            $(event.delegateTarget).addClass("selected");
            $this.selected = true;
            $this.value = $this.selected
        }
        
        if ($this.callBack != undefined ) {
            $this.callBack($this.callBackParams);
        }
    }
});

var UIText = UIEle.$extend({
     __init__ : function(dom, callBack, callBackParams, enabled) {
        this.$super(dom, enabled);
        
        this.name = this.dom_jq.attr("name");
        this.label_text = this.dom_jq.siblings("label").text();
        this.value = undefined;
        this.default_value = this.dom_jq.val();
        
        this.callBack = callBack;
        this.callBackParams = callBackParams;
        ////console.debug(this.dom_jq);
        
        this.dom_jq.keyup({$this:this}, this.textChange);
        //this.dom_jq.focus({$this:this}, this.textFocus);
        this.dom_jq.click({$this:this}, this.textFocus);
        
     },
     
     textChange : function(event) {
        $this = event.data.$this
        $this.value = $this.dom_jq.val();
        
        if ( $this.callBack ) {
            $this.callBack($this.callBackParams);
        }
     },
     
     textFocus : function(event) {
        $this = event.data.$this
        
        if ( $this.dom_jq.val() == $this.default_value ) {
            $this.dom_jq.select();
        }
    }
     
});


var UIModalDialog = UIEle.$extend({
    __init__ : function(dom, server_data) {
        this.$super(dom);
        
        this.form_elements = [];
        
        $(window).resize({$this:this}, this.resizeWindow);
        
        this.fade_duration = 100;
        
        this.server_data = server_data;
        
        this.showDialog();
        
        // need to find all of the form elements in the active dialog
        
        this.action = this.active_dialog.find("form").attr("action");
        
        // we always have these two buttons
        this.go_button = new UIButton(this.active_dialog.find( "[name='go']" ), this.doGo, {$this:this}, false);
        this.form_elements.push( new UIButton(this.active_dialog.find( "[name='cancel']" ), this.doCancel, {$this:this}, true) );
        this.form_elements.push( this.go_button );
        
        $this = this;
        //look for checkboxes
        this.active_dialog.find( "[type='checkbox']" ).each( function (index) {
            $this.form_elements.push( new UICheckbox(this, $this.activitySelected, {$this:$this}) );
        });
        
        //look for textfields
        this.active_dialog.find( "[type='text']" ).each( function (index) {
            $this.form_elements.push( new UIText(this) );
        });
        
        $(window).keyup({$this:this}, this.checkReturn);
        
    },
    
    checkReturn : function(event) {
        $this = event.data.$this;
        
        if ( event.which == 13 ) {
            if( $this.go_button.enabled ) {
                $this.doGo({$this:$this});
            }
        }
    },
    
    activitySelected : function(params) {
        $this = params.$this;
        
        
        activate_go = false;
        for ( element in $this.form_elements ) {
            if ( $this.form_elements[element] instanceof UICheckbox ) {
                if ( $this.form_elements[element].value ) {
                    activate_go = true;
                }
            }
        }
        ////console.debug(activate_go);
        if ( activate_go ) {
            $this.go_button.UIEleEnable($this.go_button);
        } else {
            $this.go_button.UIEleDisable($this.go_button);
        }
    },
    
    doCancel : function(params) {
        $this = params.$this;
        
        $this.hideDialog();
    },
    
    doGo : function(params) {
        $this = params.$this;
        
        query_str = "?";
        query_obj = {};
        for ( element in $this.form_elements ) {
            if ( $this.form_elements[element] instanceof UICheckbox ) {
                query_str += $this.form_elements[element].name + "=" + $this.form_elements[element].value + "&";
                query_obj[$this.form_elements[element].name] = $this.form_elements[element].value
            }
            if ( $this.form_elements[element] instanceof UIText ) {
                if ( $this.form_elements[element].value != undefined ) {
                    query_str += $this.form_elements[element].name + "=" + $this.form_elements[element].value + "&";
                    query_obj[$this.form_elements[element].name] = $this.form_elements[element].value
                }
            }
        }
        query_obj.server_data = JSON.stringify($this.server_data);
        ////console.debug(JSON.stringify($this.server_data));
        //query_str = query_str.substr(0, query_str.length - 1);
        
        /* this is an overlay but for some reason it doesnt work.
        $this.activity_overlay = $("<div class=\"activity_overlay\"><div class=\"activity_overlay_message\">Please wait.</div></div>").appendTo("#uber");
        $this.activity_overlay.css({"position":"absolute", "display":"none"});
        $this.activity_overlay.offset(
            {   top:    0,
                left:   0
        });
        $this.activity_overlay.width($(window).width());
        $this.activity_overlay.height($("html").height());
        $this.activity_overlay.fadeIn($this.fade_duration, function() {
            document.location.href = $this.action + query_str + "server_data=" + JSON.stringify($this.server_data);
        });
        */
        
        //console.debug($this.action + query_str + "server_data=" + JSON.stringify($this.server_data))
        
        $this.activity_overlay = $("<div class=\"activity_overlay\"><div class=\"activity_overlay_message\">Please wait.</div></div>").appendTo("#uber");
        $this.activity_overlay.css({"position":"absolute", "display":"none"});
        $this.activity_overlay.offset(
            {   top:    0,
                left:   0
        });
        $this.activity_overlay.width($(window).width());
        $this.activity_overlay.height($("html").height());
        $this.activity_overlay.fadeIn($this.fade_duration);
        
        $.get($this.action, query_obj, function(data) {
            ////console.debug([$this, $this.activity_overlay]);
            //$this.activity_overlay.remove();
            newDoc = document.open("text/html");
            newDoc.write(data);
            newDoc.close();
            //$(window).html(data);
            //$this.detailLoadComplete($this);
        });
        
        //document.location.href = $this.action + query_str + "server_data=" + JSON.stringify($this.server_data);
        
        ////console.debug([$this.action, query_str]);
    },
    
    showDialog : function() {
        // create an overlay to capture all events outside the dialog box
        this.overlay = $("<div class=\"dialog_overlay\"></div>").appendTo("#uber");
        
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
        $(window).off("keyup");
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
        
        // if there is an activity overlay, it needs to resize
        if ( $this.activity_overlay != undefined ) {
            $this.activity_overlay.width(the_width);
            $this.activity_overlay.height($("html").height());
        }
    }

});

var UISingleSelectPopup = UIPopup.$extend({
    __init__ : function(dom, callBack, callBackParams, active) {
        this.$super(dom, active);
        
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