var VCAction = Class.$extend({

__init__ : function(dom_id, back_action) {
    this.dom_id = dom_id;
    this.dom_jq = $(this.dom_id);
    
    this.back_action = back_action;
    
    this.dom_jq.prepend("<div id=\"back\"></div>");
    this.back_button = $("#back");
    
    this.back_button.mouseenter({$this:this}, this.backHovered);
    this.back_button.mouseleave({$this:this}, this.backUnHovered);
    this.back_button.click({$this:this}, this.backClicked);
    
    this.resizeWindow();
    $(window).resize({$this:this}, this.resizeWindow);
},

resizeWindow : function(event) {
    // the FaultMap object
    if ( event == undefined ) {
        $this = this
    } else {
        $this = event.data.$this
    }
    
    the_width = $(window).outerWidth();
    
    content_height = $(window).height() - $("#system").outerHeight(true) - $("#menu_bar").outerHeight(true) - parseFloat($("#uber").css("padding-top").split("px")[0]);
    
    $this.dom_jq.height(content_height);
    
    top_items_height = Math.max($("#action_info").outerHeight(true) + $("#action_time_series").outerHeight(true) + $this.back_button.outerHeight(true), $("#action_fault_map").outerHeight(true) + $this.back_button.outerHeight(true));
    
    $("#action_fault_map").css({"position":"absolute", "top":$this.back_button.outerHeight(true), "left":$("#action_time_series").outerWidth(true)});
    
    ////console.debug(top_items_height);
    
    scaling_row = 0;
    $(".action_scaling_plot").each(function(index) {
        
        ////console.debug(parseFloat($(this).css("margin-top").split("px")[0]) + parseFloat($this.back_button.css("margin-bottom").split("px")[0]));
        ////console.debug($(this).outerHeight(true) - $(this).outerHeight());
        
        $(this).css({"position":"absolute",
            "top":top_items_height + scaling_row * $(this).outerHeight(true),
            "left":(index%2) * ($("#action_time_series").outerWidth() + parseFloat($("#action_fault_map").css("margin-left").split("px")[0]))
        });
        if ( index%2 == 1 )
            scaling_row += 1;
    });
    
    //if (scaling_row == 0)
    //    scaling_row = 1;
    
    console.debug(scaling_row);
    st_row = 0;
    $(".action_space_time_plot").each(function(index) {
    
        $(this).css({"position":"absolute",
            "top":top_items_height + scaling_row * $(".action_scaling_plot").outerHeight(true) + st_row * $(this).outerHeight(true)
        });
        if ( index%2 == 1 )
            st_row += 1;
    });

    
},

backClicked : function(event) {
    // the EventList object
    $this = event.data.$this
    
    document.location.href = $this.back_action;
    
    //if ( $this.from_view == "fault_map" )
    //    document.location.href =
    ////console.debug($this.from_view);
},

backHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    //if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).addClass("hovered");
    //}
    $(event.delegateTarget).css("cursor","pointer");
},

backUnHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    //if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).removeClass("hovered");
    //}
    $(event.delegateTarget).css("cursor","auto");
}

});