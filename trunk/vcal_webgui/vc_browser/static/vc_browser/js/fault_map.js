var FaultMap = Class.$extend({

__init__ : function(map_dom_id, list_dom_id, fault_map_data) {
    this.map_dom_id = map_dom_id;
    this.map_dom_jq = $(this.map_dom_id);
    
    this.list_dom_id = list_dom_id;
    this.list_dom_jq = $(this.list_dom_id);
    
    this.fault_map_data = jQuery.parseJSON(fault_map_data);
    
    ////console.debug(this.fault_map_data);
    
    //this.sections = sections;
    this.open_info_window = {};
    //this.selected_sections = {};
    this.number_of_sections = 0;
    
    this.server_data = {from_view: 'fault_map', selected:{}, unselected:{}};
    
    this.selected_line_inner_color = "#FFFFFF";
    this.unselected_line_inner_color = "#FF0000";
    
    this.unselected_line_options = {
        strokeColor: this.unselected_line_inner_color,
        strokeOpacity: 1.0,
        strokeWeight: 5
    };
    this.selected_line_options = {
        strokeColor: "#000000",
        strokeOpacity: 1.0,
        strokeWeight: 10
    };
    
    this.list_button_template = "<div id=\"fault_list_button\"></div>";
    this.list_item_template = "<li>{{sid}} {{sname}}</li>";
    this.infoWindow_template = "<div class=\"info_window\">{{sid}} {{sname}}</div>";
    this.list_visible = false;
    this.slide_duration = 100;
    
    // need this to access the object from within jquery loops
    $this = this;
    
    // initialize the fault list button
    this.list_dom_jq.after(this.list_button_template);
    this.list_button_dom_jq = $("#fault_list_button");
    lb_top = this.list_dom_jq.height() - 6;
    this.list_button_dom_jq.css({"position":"absolute", "left":"0px"});
    //this.list_button_dom_jq.mouseenter({$this:this}, this.listButtonHovered);
    //this.list_button_dom_jq.mouseleave({$this:this}, this.listButtonUnHovered);
    this.list_button_dom_jq.click({$this:this}, this.listButtonClicked);
    
    // initialize the fault list controls
    new UIButton($("#fault_list_header").children( "[name='all']" ), this.faultMapSelectAllCallback, {$this:this});
    new UIButton($("#fault_list_header").children( "[name='none']" ), this.faultMapSelectNoneCallback, {$this:this});
    
    this.search_field = new UIText(this.list_dom_jq.find("input[name=\"search\"]"), this.searchFieldChange, {$this:this});
    
    //this.all_button = this.list_dom_jq.find("input[name=\"all\"]");
    //this.none_button = this.list_dom_jq.find("input[name=\"none\"]");
    //this.search_field = this.list_dom_jq.find("input[name=\"search\"]");
    
    //this.all_button.mouseenter({$this:this}, this.allButtonHovered);
    //this.all_button.mouseleave({$this:this}, this.allButtonUnHovered);
    //this.all_button.click({$this:this}, this.allButtonClicked);
    
    //this.none_button.mouseenter({$this:this}, this.noneButtonHovered);
    //this.none_button.mouseleave({$this:this}, this.noneButtonUnHovered);
    //this.none_button.click({$this:this}, this.noneButtonClicked);
    
    //this.search_field.keyup({$this:this}, this.searchFieldChange);
    
    //set the height of the map window
    this.resizeWindow();
    $(window).resize({$this:this}, this.resizeWindow);
    
    // initilize the map and the list
    mapOptions = {
        center: new google.maps.LatLng(this.fault_map_data.system_center[0], this.fault_map_data.system_center[1]),
        mapTypeId: google.maps.MapTypeId.TERRAIN,
        disableDefaultUI: true,
        scaleControl: true,
        zoomControl: true,
        zoomControlOptions: {
            position: google.maps.ControlPosition.TOP_RIGHT
        },
        mapTypeControl: true
    };
    
    this.map = new google.maps.Map(this.map_dom_jq[0], mapOptions);
    
    this.list_dom_jq.children("#fault_list_list").append("<ul></ul>");
    
    all_sections = new google.maps.LatLngBounds( );
    
    //google.maps.event.addListener(this.map, "click", this.mapClick);
    
    this.sections = {};
    sids = [];
    
    for (var key in this.fault_map_data) {
        if (this.fault_map_data.hasOwnProperty(key)) {
            if ( key != 'system_center' ) {
                this.sections[key] = this.fault_map_data[key]
                sids.push(parseInt(key));
            }
        }
    }
    
    sids.sort(function(a,b){return a-b});
    
    len = sids.length;
    
    for (i = 0; i < len; i++ ) {
        
        sid = "" + sids[i];
        this_section = [];
        for ( j = 0; j < this.sections[sid].trace.length; j++) {
            this_section.push(new google.maps.LatLng(this.sections[sid].trace[j][0], this.sections[sid].trace[j][1]));
            all_sections.extend(new google.maps.LatLng(this.sections[sid].trace[j][0], this.sections[sid].trace[j][1]));
        }
        
        //set up the unselected line and add it to the map
        this.unselected_line_options.path = this_section
        this.unselected_line_options.zIndex = sid;
        
        this.sections[sid].unselected_polyLine = new google.maps.Polyline(this.unselected_line_options);
        this.sections[sid].unselected_polyLine.section_name = this.sections[sid].name;
        this.sections[sid].unselected_polyLine.section_id = sid;
        this.sections[sid].unselected_polyLine.$this = $this;
        this.sections[sid].unselected_polyLine.setMap(this.map);
        google.maps.event.addListener(this.sections[sid].unselected_polyLine, 'click', this.sectionClick);
        google.maps.event.addListener(this.sections[sid].unselected_polyLine, 'mouseover', this.sectionHover);
        google.maps.event.addListener(this.sections[sid].unselected_polyLine, 'mouseout', this.sectionUnHover);
        
        //set up the selected line but dont add it to the map yet
        this.selected_line_options.path = this_section
        
        this.sections[sid].selected_polyLine = new google.maps.Polyline(this.selected_line_options);
        this.sections[sid].selected_polyLine.section_name = this.sections[sid].name;
        this.sections[sid].selected_polyLine.section_id = sid;
        this.sections[sid].selected_polyLine.$this = $this;
        google.maps.event.addListener(this.sections[sid].selected_polyLine, 'click', this.sectionClick);
        google.maps.event.addListener(this.sections[sid].selected_polyLine, 'mouseover', this.sectionHover);
        google.maps.event.addListener(this.sections[sid].selected_polyLine, 'mouseout', this.sectionUnHover);
        
        //set up the info window
        content_string = sid + " " + this.sections[sid].name;
        this.sections[sid].infoWindow = new google.maps.InfoWindow({
            content: this.infoWindow_template.replace("{{sid}}", sid).replace("{{sname}}", this.sections[sid].name),
            position:this_section[0],
            disableAutoPan:true,
            maxWidth:100
        });
        
        //add this section to the fault_list
        this.list_dom_jq.children("#fault_list_list").children("ul").append(this.list_item_template.replace("{{sid}}", sid).replace("{{sname}}", this.sections[sid].name));
        
        ////console.debug(this.list_dom_jq.children("#fault_list_list").children("ul").children("li").last());
        
        this.sections[sid].list_item = this.list_dom_jq.children("#fault_list_list").children("ul").children("li").last();
        this.sections[sid].list_item.data("section_id", sid);
        
        // set the list item call back functions
        this.sections[sid].list_item.mouseenter({$this:this, section_id:sid}, this.listItemHovered);
        this.sections[sid].list_item.mouseleave({$this:this, section_id:sid}, this.listItemUnHovered);
        this.sections[sid].list_item.click({$this:this, section_id:sid}, this.listItemClicked);
        
        this.number_of_sections += 1;
            
    }
    
    this.map.fitBounds( all_sections );
    
},

resizeWindow : function(event) {
    // the FaultMap object
    if ( event == undefined ) {
        $this = this
    } else {
        $this = event.data.$this
    }
    
    content_area_height = $(window).height() - $("#system").outerHeight(true) - $("#menu_bar").outerHeight(true) - parseFloat($("#uber").css("padding-top").split("px")[0]);
    //content_area_height = $(window).height() - 41 - 28 - 10;
    //list_height = content_area_height  - 35;
    list_height = content_area_height  - $this.list_dom_jq.children("#fault_list_header").outerHeight(true);
    //list_height = content_area_height  - $this.list_dom_jq.children("#fault_list_header").outerHeight(true);// - parseFloat($this.list_dom_jq.children("#fault_list_list").css("padding-top").split("px")[0]) - parseFloat($this.list_dom_jq.children("#fault_list_list").css("padding-bottom").split("px")[0]);
    
    ////console.debug(parseFloat($this.list_dom_jq.children("#fault_list_list").css("padding-top").split("px")[0]));
    $this.map_dom_jq.height(content_area_height);
    $this.list_dom_jq.height(content_area_height);
    $this.list_dom_jq.children("#fault_list_list").height(list_height);
    $this.list_button_dom_jq.css({"top":(content_area_height/2.0 - $this.list_button_dom_jq.height()/2.0) + "px"});

},

searchFieldChange : function(params) {
    // the FaultMap object
    $this = params.$this
    
    curr_search_val = new RegExp( $this.search_field.value, "i" )
    
    //console.debug($this.search_field.value);
    
    for (var sid in $this.sections) {
        if ( $this.sections[sid].name.search(curr_search_val) == 0 ) {
            $this.sections[sid].list_item.css({"display":"block"});
            ////console.debug($this.sections[sid].name);
        } else {
            $this.sections[sid].list_item.css({"display":"none"});
        }
        //$this.sections[sid].name;
    }
    
},

allButtonHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    
    //if ( !$this.all_selected ) {
    //this item has not been selected
    $(this).addClass("hovered");
    //$(this).css({"background-color":"#FF0000", "color":"#FFFFFF"});
    //}
    $(event.delegateTarget).css("cursor","pointer");
},

allButtonUnHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    
    //if ( !$this.all_selected ) {
    //this item has not been selected
    $(this).removeClass("hovered");
    //$(this).css({"background-color":"#FFFFFF", "color":"#000000"});
    //}
    $(event.delegateTarget).css("cursor","auto");
},

faultMapSelectAllCallback : function(params) {
    // the FaultMap object
    $this = params.$this;
    
    //if ( !$this.all_selected ) {
    // select everything
    for (var sid in $this.sections) {
        if ( $this.sections[sid].list_item.is(":visible") )
            $this.server_data.selected[sid] = true;
    }
        //$this.all_selected = true;
    //} else {
        // deselect everything
        //for (var sid in $this.sections) {
            //$this.server_data.selected[sid] = false;
        //}
        //$this.all_selected = false;
    //}
    
    $this.setSelectedSectionsProperties();
},

noneButtonHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    
    //if ( !$this.all_selected ) {
    //this item has not been selected
    $(this).addClass("hovered");
    //$(this).css({"background-color":"#FF0000", "color":"#FFFFFF"});
    //}
    $(event.delegateTarget).css("cursor","pointer");
},

noneButtonUnHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    
    //if ( !$this.all_selected ) {
    //this item has not been selected
    $(this).removeClass("hovered");
    //$(this).css({"background-color":"#FFFFFF", "color":"#000000"});
    //}
    $(event.delegateTarget).css("cursor","auto");
},

faultMapSelectNoneCallback : function(params) {
    // the FaultMap object
    $this = params.$this;
    
    //if ( !$this.all_selected ) {
    // select everything
    //for (var sid in $this.sections) {
    $this.server_data.selected = {};
    //}
        //$this.all_selected = true;
    //} else {
        // deselect everything
        //for (var sid in $this.sections) {
            //$this.server_data.selected[sid] = false;
        //}
        //$this.all_selected = false;
    //}
    
    $this.setSelectedSectionsProperties();
},


listItemHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    section_id = event.data.section_id;
    
    if ( $this.server_data.selected[section_id] == undefined || !$this.server_data.selected[section_id] ) {
        //this item has not been selected
        $(this).addClass("hovered");
        //$(this)
        //$(this).css({"background-color":"#FF0000", "color":"#FFFFFF"});
    }
    $(event.delegateTarget).css("cursor","pointer");
},

listItemUnHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    section_id = event.data.section_id;
    
    if ( $this.server_data.selected[section_id] == undefined || !$this.server_data.selected[section_id] ) {
        //this item has not been selected
        $(this).removeClass("hovered");
        //$(this).css({"background-color":"#FFFFFF", "color":"#000000"});
    }
    $(event.delegateTarget).css("cursor","auto");
},

listItemClicked : function(event) {
    // the FaultMap object
    $this = event.data.$this
    section_id = event.data.section_id;
    
    ////console.debug(event);
    
    if ( event.metaKey ) {
        if ( $this.server_data.selected[section_id] == undefined || !$this.server_data.selected[section_id] ) {
            $this.server_data.selected[section_id] = true;
        }
        else {
            $this.server_data.selected[section_id] = false;
        }
    }
    else {
        
        for (var sid in $this.server_data.selected) {
            if ($this.server_data.selected.hasOwnProperty(sid)) {
                $this.server_data.selected[sid] = false;
            }
        }
        
        $this.server_data.selected[section_id] = true;
    }
    
    
    $this.setSelectedSectionsProperties();
    
    //$this.server_data.selected[this.section_id] = true;
    //$this.setSelectedSectionsProperties();
    //$(this).css({"background-color":"#808080", "color":"#000000"});
},

listButtonHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    
    if ($this.list_visible) {
        $this.list_button_dom_jq.css({"background-image":"url(images/left_arrow_active.png)"});
    }
    else {
        $this.list_button_dom_jq.css({"background-image":"url(images/right_arrow_active.png)"});
    }
},

listButtonUnHovered : function(event) {
    // the FaultMap object
    $this = event.data.$this
    
    if ($this.list_visible) {
        $this.list_button_dom_jq.css({"background-image":"url(images/left_arrow.png)"});
    }
    else {
        $this.list_button_dom_jq.css({"background-image":"url(images/right_arrow.png)"});
    }
},

listButtonClicked : function(event) {
    // the FaultMap object
    $this = event.data.$this
    
    if ($this.list_visible) {
        //$this.list_button_dom_jq.css({"background-image":"url(images/right_arrow.png)"});
        $this.list_button_dom_jq.animate({
            left:"0px"
        }, $this.slide_duration
        );
        $this.list_dom_jq.animate({
            left: "-" + ($this.list_dom_jq.width() + 1) + "px"
        }, $this.slide_duration, function() {
            $this.list_visible = false;
            $(this).hide();
        });
        //$this.list_dom_jq.hide();
    }
    else {
        //$this.list_button_dom_jq.css({"background-image":"url(images/left_arrow.png)"});
        $this.list_dom_jq.show();
        $this.list_button_dom_jq.animate({
            left:($this.list_dom_jq.width()) + "px"
        }, $this.slide_duration
        );
        $this.list_dom_jq.animate({
            left: "0px"
        }, $this.slide_duration, function() {
            $this.list_visible = true;
        });
    }
},

mapClick : function(event) {
    // the FaultMap object
    $this = this.$this
    
    //console.debug(this);
},

sectionHover : function(event) {
    // the FaultMap object
    $this = this.$this
    
    $this.sections[this.section_id].infoWindow.open($this.map);
    
    $this.sections[this.section_id].unselected_polyLine.setOptions({ zIndex: $this.number_of_sections + 2 });
    $this.sections[this.section_id].selected_polyLine.setOptions({zIndex: $this.number_of_sections + 1 });
    $this.sections[this.section_id].selected_polyLine.setMap($this.map);
},

sectionUnHover : function(event) {
    // the FaultMap object
    $this = this.$this
    
    $this.sections[this.section_id].infoWindow.close();
    
    if ( $this.server_data.selected[this.section_id] == undefined || !$this.server_data.selected[this.section_id] ) {
        $this.sections[this.section_id].selected_polyLine.setMap(null);
        $this.sections[this.section_id].unselected_polyLine.setOptions({ zIndex: this.section_id });
    }
},

sectionClick : function(event) {
    
    // the FaultMap object
    $this = this.$this
    
    function type(o){
        return !!o && Object.prototype.toString.call(o).match(/(\w+)\]/)[1];
    }
    
    //for ( prop in event ) {
    //    //console.debug([prop, event[prop]]);
    //}
    
    for ( level1 in event ) {
        ////console.debug( type( event[level1] ) );
        if ( type( event[level1] ) == "MouseEvent" ) {
            MouseEvent = event[level1];
            break;
        }
        MouseEvent = undefined;
    }
    ////console.debug(event);
    ////console.debug(findObjectByLabel(event, "metaKey"));
    
    ////console.debug(MouseEvent);
    
    if ( MouseEvent.metaKey ) {
        if ( $this.server_data.selected[this.section_id] == undefined || !$this.server_data.selected[this.section_id] ) {
            $this.server_data.selected[this.section_id] = true;
        }
        else {
            $this.server_data.selected[this.section_id] = false;
        }
    }
    else {
        
        for (var sid in $this.server_data.selected) {
            if ($this.server_data.selected.hasOwnProperty(sid)) {
                $this.server_data.selected[sid] = false;
            }
        }
        
        $this.server_data.selected[this.section_id] = true;
    }
    
    $this.setSelectedSectionsProperties();
},

setSelectedSectionsProperties : function() {
    if ( jQuery.isEmptyObject( this.server_data.selected ) ) {
        //nothing is selected: unselect eveything
        for (var sid in this.sections) {
            if (this.sections.hasOwnProperty(sid)) {
                this.sections[sid].selected_polyLine.setMap(null);
                    
                this.sections[sid].unselected_polyLine.setOptions({
                    zIndex:sid,
                    strokeColor:this.unselected_line_inner_color
                });
                this.sections[sid].list_item.removeClass("selected hovered");
            }
        }
    } else {
        // something is selected
        for (var sid in this.server_data.selected) {
            if (this.server_data.selected.hasOwnProperty(sid)) {
                ////console.debug(sid + " " + this.server_data.selected[sid] + " " + this.sections[sid].name);
                if (this.server_data.selected[sid]) {
                    this.sections[sid].unselected_polyLine.setOptions({
                        zIndex:this.number_of_sections + 2,
                        strokeColor:this.selected_line_inner_color
                    });
                    this.sections[sid].selected_polyLine.setOptions({zIndex: this.number_of_sections + 1});
                    this.sections[sid].selected_polyLine.setMap(this.map);
                    this.sections[sid].list_item.addClass("selected");
                    //this.sections[sid].list_item.css({"background-color":"#808080", "color":"#000000"});
                }
                else {
                    this.sections[sid].selected_polyLine.setMap(null);
                    
                    this.sections[sid].unselected_polyLine.setOptions({
                        zIndex:sid,
                        strokeColor:this.unselected_line_inner_color
                    });
                    this.sections[sid].list_item.removeClass("selected hovered");
                    //this.sections[sid].list_item.css({"background-color":"#FFFFFF", "color":"#000000"});
                }
            }
        }
    }
}

});