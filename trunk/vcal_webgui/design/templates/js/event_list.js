var EventList = Class.$extend({

__init__ : function(_el_dom_id, _ed_dom_id, initial_list) {
    $this = this;
    
    this._el_dom_id = _el_dom_id;
    this._el_dom_jq = $(this._el_dom_id);
    
    this._ed_dom_id = _ed_dom_id;
    this._ed_dom_jq = $(this._ed_dom_id);
    
    this.server_data = {selected:{}, unselected:{}};
    
    this.slide_duration = 100;
    
    this.list_dat = jQuery.parseJSON(initial_list);
    
    this._el_table_dom_jq = null;
    
    this.number_items_showing = 500;
    
    this.current_page = 1;
    this.first_page = 1;

    
    
    this.max_pages = 40;
    
    this.column_info = {
        evid:{
            column_name:"Number",
            column_width:8,
            sort:"ascending"
            },
        evyear:{
            column_name:"Year",
            column_width:15
            },
        evmag:{
            column_name:"Magnitude",
            column_width:15
            },
        evtriggersec:{
            column_name:"Trigger Section",
            column_width:14
            },
        evtriggerele:{
            column_name:"Trigger Element",
            column_width:10
            },
        evinvolvedsec:{
            column_name:"# of Sections Involved",
            column_width:10
            },
        evinvolvedele:{
            column_name:"# of Elements Involved",
            column_width:10
            },
        evaveslip:{
            column_name:"Average Slip [m]",
            column_width:15
            },
        evdetail:{
            column_name:"",
            column_width:3
            }
    };
    
    this._el_dom_jq.append("<div id=\"event_list_controls\"></div>");
    this._el_controls_dom_jq = this._el_dom_jq.children("#event_list_controls");
    
    this._el_controls_dom_jq.append("<div id=\"event_list_selection\">\
                                        <input class=\"button small\" name=\"all\" type=\"button\" value=\"All\" />\
                                        <input class=\"button small\" name=\"none\" type=\"button\" value=\"None\" />\
                                    </div>\
                                    <div id=\"event_list_pager\"></div>\
                                    <div id=\"event_list_display_items\" class=\"ss_popup\">\
                                        <ul>\
                                            <li>10</li>\
                                            <li>50</li>\
                                            <li>100</li>\
                                            <li class=\"selected\">500</li>\
                                            <li>1000</li>\
                                        </ul>\
                                    </div>");
    
    // find the widths for the pager
    // first one character and the arrows
    $("#event_list_pager").append("<div class=\"page\">0</div>");
    $("#event_list_pager").append("<div class=\"arrow_left\"></div>");
    
    this.pageWidth_1char = $(".page").outerWidth();
    this.pagerArrowWidth = $(".arrow_left").outerWidth();
    
    $(".page, .arrow_left").remove();
    
    // now with 2 chars
    $("#event_list_pager").append("<div class=\"page\">00</div>");
    this.pageWidth_2char = $(".page").outerWidth();
    $(".page").remove();
    
    // now with 3 chars
    $("#event_list_pager").append("<div class=\"page\">000</div>");
    this.pageWidth_3char = $(".page").outerWidth();
    $(".page").remove();
    
    // now with 4 chars
    $("#event_list_pager").append("<div class=\"page\">0000</div>");
    this.pageWidth_4char = $(".page").outerWidth();
    $(".page").remove();
    
    // now with 5 chars. I hope this is all we need
    $("#event_list_pager").append("<div class=\"page\">00000</div>");
    this.pageWidth_5char = $(".page").outerWidth();
    $(".page").remove();
    
    this.setPageWidth();
    
    new UISingleSelectPopup("#event_list_display_items", this.eventListDisplayItemCallback, {$this:this});
    new UIButton($("#event_list_selection").children( "[name='all']" ), this.eventListSelectAllCallback, {$this:this});
    new UIButton($("#event_list_selection").children( "[name='none']" ), this.eventListSelectNoneCallback, {$this:this});
    
    this.buildList();
    this.resizeWindow();
    $(window).resize({$this:this}, this.resizeWindow);
},

// this needs to be called every time the number of items showing is updated 
setPageWidth : function() {
    this.first_page = 1;
    this.total_pages = Math.ceil(this.list_dat.TotalRecordCount/this.number_items_showing);
    //this.total_pages = 1000;
    length = (this.total_pages + "").length;
    switch (length) {
        case 1:
            this.pageWidth = this.pageWidth_1char;
            break;
        case 2:
            this.pageWidth = this.pageWidth_2char;
            break;
        case 3:
            this.pageWidth = this.pageWidth_3char;
            break;
        case 4:
            this.pageWidth = this.pageWidth_4char;
            break;
        case 5:
            this.pageWidth = this.pageWidth_5char;
            break;
    }
},

selectAll : function () {
    this.server_data.selected = {all : true};
    this.server_data.unselected = {};
    $("#event_list tr").addClass("selected");
    //console.debug([this.server_data.selected, this.server_data.unselected])
},

selectNone : function () {
    this.server_data.selected = {};
    this.server_data.unselected = {all : true};
    $("#event_list tr").removeClass("selected");
    //console.debug([this.server_data.selected, this.server_data.unselected])
},

eventListSelectAllCallback : function(params) {
    $this = params.$this;
    
    $this.selectAll();
},

eventListSelectNoneCallback : function(params) {
    $this = params.$this;
    
    $this.selectNone();
},

buildPager : function() {
    pager_jq = $("#event_list_pager");
    
    pager_jq.empty();
    
    this.total_pages = Math.ceil(this.list_dat.TotalRecordCount/this.number_items_showing);
    
    if ( this.total_pages > this.max_pages ) {
        showing_pages = this.max_pages;
    } else {
        showing_pages = this.total_pages;
    }
    
    
    // create the pager content
    pager_content = "";
    
    // open the pager table. Yes we are using a table. Get over it.
    pager_content += "<table class=\"pages\"><tr>"
    
    // left arrows
    pager_content += "<td class=\"double_arrow_left\"></td><td class=\"arrow_left\"></td>";
    
    // The pages
    for ( i = this.first_page; i < this.first_page + this.max_pages; i++ ) {
        pager_content += "<td class=\"page\" style=\"width:" + this.pageWidth + "px;\">" + i + "</td>";
    }
    
    // The ... cell at the end to indicate there are more pages to the right.
    if ( this.first_page + this.max_pages < this.total_pages ) {
        pager_content += "<td class=\"page\" style=\"width:" + this.pageWidth + "px;\">...</td>";
    } else {
        pager_content += "<td class=\"page\" style=\"width:" + this.pageWidth + "px;\">   </td>";
    }
    
    // right arrows
    pager_content += "<td class=\"arrow_right\"></td><td class=\"double_arrow_right\"></td>";
    
    // close the table
    pager_content += "</tr></table>";
    
    // add the table to the dom
    pager_jq.append(pager_content);
    
    // decite whether to show the arrows
    if ( this.first_page > this.max_pages ) {
        $("#event_list_pager .double_arrow_left").addClass(
            "active"
        ).mouseenter(
            {$this:this},
            $this.pagerArrowHovered
        ).mouseleave(
            {$this:this},
            $this.pagerArrowUnHovered
        ).click(
            {$this:this, arrow:"double_left"},
            $this.pagerArrowClicked
        );
        $("#event_list_pager .arrow_left").addClass(
            "active"
        ).mouseenter(
            {$this:this},
            $this.pagerArrowHovered
        ).mouseleave(
            {$this:this},
            $this.pagerArrowUnHovered
        ).click(
            {$this:this, arrow:"left"},
            $this.pagerArrowClicked
        );

    } else {
        $("#event_list_pager .double_arrow_left, #event_list_pager .arrow_left").removeClass(
            "active"
        ).off(
            "mouseenter"
        ).off(
            "mouseleave"
        ).off(
            "click"
        );
    }
    
    if ( this.first_page + this.max_pages < this.total_pages ) {
        $("#event_list_pager .arrow_right").addClass(
            "active"
        ).mouseenter(
            {$this:this},
            $this.pagerArrowHovered
        ).mouseleave(
            {$this:this},
            $this.pagerArrowUnHovered
        ).click(
            {$this:this, arrow:"right"},
            $this.pagerArrowClicked
        );
        $("#event_list_pager .double_arrow_right").addClass(
            "active"
        ).mouseenter(
            {$this:this},
            $this.pagerArrowHovered
        ).mouseleave(
            {$this:this},
            $this.pagerArrowUnHovered
        ).click(
            {$this:this, arrow:"double_right"},
            $this.pagerArrowClicked
        );
    } else {
        $("#event_list_pager .double_arrow_right, #event_list_pager .arrow_right").removeClass(
            "active"
        ).off(
            "mouseenter"
        ).off(
            "mouseleave"
        ).off(
            "click"
        );
    }

    // callback functions for the pages
    $("#event_list_pager .pages .page").mouseenter(
        {$this:this},
        $this.pageHovered
    ).mouseleave(
        {$this:this},
        $this.pageUnHovered
    ).click(
        {$this:this},
        $this.pageClicked
    );

},

pagerArrowClicked : function(event) {
    // the EventList object
    $this = event.data.$this;
    arrow = event.data.arrow;
    
    length = (this.total_pages + "").length;
    
    //console.debug([$this.max_pages, $this.total_pages])
    
    switch ( arrow ) {
        case "left":
            $this.first_page -= $this.max_pages;
            break;
        case "double_left":
            $this.first_page -= $this.max_pages * length;
            break;
        case "right":
            $this.first_page += $this.max_pages;
            break;
        case "double_right":
            $this.first_page += $this.max_pages * length;
            break;
    }
    
    if ( $this.first_page < 1 )
        $this.first_page = 1;
        
    if ( $this.first_page + $this.max_pages > $this.total_pages )
        $this.first_page = $this.total_pages - $this.max_pages + 1;
    
    $this.buildPager();
    $this.resizeWindow();
    
},

pageClicked : function(event) {
    // the EventList object
    $this = event.data.$this;
    
    selected_page = parseInt($(event.delegateTarget).text());
    
    //console.debug(selected_page);
},

pageUnHovered : function(event) {
    // the EventList object
    //$this = event.data.$this
    
    $(event.delegateTarget).removeClass("hovered");
    $(event.delegateTarget).css("cursor","auto");
},

pageHovered : function(event) {
    // the EventList object
    //$this = event.data.$this
    
    $(event.delegateTarget).addClass("hovered");
    $(event.delegateTarget).css("cursor","pointer");
},

pagerArrowUnHovered : function(event) {
    // the EventList object
    //$this = event.data.$this
    
    $(event.delegateTarget).removeClass("hovered");
    $(event.delegateTarget).css("cursor","auto");
},

pagerArrowHovered : function(event) {
    // the EventList object
    //$this = event.data.$this
    
    $(event.delegateTarget).addClass("hovered");
    $(event.delegateTarget).css("cursor","pointer");
},


eventListDisplayItemCallback : function(params) {
    $this = params.$this;
    
    //console.debug("event_list_test_"+params.selection+"_dat.json");
    
    $this.number_items_showing = parseInt(params.selection);
    
    for ( column in $this.column_info ) {
        if ( $this.column_info[column].sort != undefined ) {
            current_sort_direction = $this.column_info[column].sort;
            current_sort_column = column;
        }
    }
    
    //console.debug("event_list_test_"+params.selection+"_dat_" + current_sort_column + "_" + current_sort_direction + ".json");

    $.getJSON("event_list_test_"+params.selection+"_dat_" + current_sort_column + "_" + current_sort_direction + ".json", function(data) {
        $this.list_dat = data;
        $this.setPageWidth();
        $this.buildList();
        $this.resizeWindow();
    });

},

buildList : function(event) {
    // the FaultMap object
    if ( event == undefined ) {
        $this = this
    } else {
        $this = event.data.$this
    }
    
    if ( $this._el_table_dom_jq != null ) {
        $this._el_table_dom_jq.remove()
    }
    
    if ( $("#event_list_header").length == 0 ) {
    
        $this._el_dom_jq.append(
            "<div id=\"event_list_header\">\
                <table border=\"0\" style=\"width:100%; padding:0px; margin:0px;\">\
                    <tr>\
                        <td class=\"event_list_header evid " + (($this.column_info.evid.sort != undefined) ? $this.column_info.evid.sort : "") + "\" style=\"width:"+ $this.column_info.evid.column_width +"%\">" + $this.column_info.evid.column_name + "</td>\
                        <td class=\"event_list_header evyear" + (($this.column_info.evyear.sort != undefined) ? $this.column_info.evyear.sort : "") + "\" style=\"width:"+ $this.column_info.evyear.column_width +"%\">" + $this.column_info.evyear.column_name + "</td>\
                        <td class=\"event_list_header evmag" + (($this.column_info.evmag.sort != undefined) ? $this.column_info.evmag.sort : "") + "\" style=\"width:"+ $this.column_info.evmag.column_width +"%\">" + $this.column_info.evmag.column_name + "</td>\
                        <td class=\"event_list_header evtriggersec" + (($this.column_info.evtriggersec.sort != undefined) ? $this.column_info.evtriggersec.sort : "") + "\" style=\"width:"+ $this.column_info.evtriggersec.column_width +"%\">" + $this.column_info.evtriggersec.column_name + "</td>\
                        <td class=\"event_list_header evtriggerele" + (($this.column_info.evtriggerele.sort != undefined) ? $this.column_info.evtriggerele.sort : "") + "\" style=\"width:"+ $this.column_info.evtriggerele.column_width +"%\">" + $this.column_info.evtriggerele.column_name + "</td>\
                        <td class=\"event_list_header evinvolvedsec" + (($this.column_info.evinvolvedsec.sort != undefined) ? $this.column_info.evinvolvedsec.sort : "") + "\" style=\"width:"+ $this.column_info.evinvolvedsec.column_width +"%\">" + $this.column_info.evinvolvedsec.column_name + "</td>\
                        <td class=\"event_list_header evinvolvedele" + (($this.column_info.evinvolvedele.sort != undefined) ? $this.column_info.evinvolvedele.sort : "") + "\" style=\"width:"+ $this.column_info.evinvolvedele.column_width +"%\">" + $this.column_info.evinvolvedele.column_name + "</td>\
                        <td class=\"event_list_header evaveslip" + (($this.column_info.evaveslip.sort != undefined) ? $this.column_info.evaveslip.sort : "") + "\" style=\"width:"+ $this.column_info.evaveslip.column_width +"%\">" + $this.column_info.evaveslip.column_name + "</td>\
                        <td class=\"event_list_header evdetail" + (($this.column_info.evdetail.sort != undefined) ? $this.column_info.evdetail.sort : "") + "\" style=\"width:"+ $this.column_info.evdetail.column_width +"%\">" + $this.column_info.evdetail.column_name + "</td>\
                    </tr>\
                </table>\
            </div>"
        );
        
        $("#event_list_header tr td").mouseenter({$this:$this}, $this.listHeaderHovered);
        $("#event_list_header tr td").mouseleave({$this:$this}, $this.listHeaderUnHovered);
        $("#event_list_header tr td").click({$this:$this}, $this.listHeaderClicked);
    }
    
    $this._el_dom_jq.append("<div id=\"event_list\"></div>");
    $this._el_table_dom_jq = $this._el_dom_jq.children("#event_list");
    
    content_str = "<table border=\"0\" style=\"width:100%\">";
    
    for ( row in $this.list_dat.Records ) {
        if ( row%2 != 0 ) {
            even_odd = "odd";
        } else {
            even_odd = "even";
        }
        
        content_str += "<tr class=\"event_list_row " + even_odd + "\">\
                            <td class=\"event_list_cell evid\" style=\"width:"+ $this.column_info.evid.column_width +"%\">" + $this.list_dat.Records[row].evid + "</td>\
                            <td class=\"event_list_cell evyear\" style=\"width:"+ $this.column_info.evyear.column_width +"%\">" + $this.list_dat.Records[row].evyear + "</td>\
                            <td class=\"event_list_cell evmag\" style=\"width:"+ $this.column_info.evmag.column_width +"%\">" + $this.list_dat.Records[row].evmag + "</td>\
                            <td class=\"event_list_cell evtriggersec\" style=\"width:"+ $this.column_info.evtriggersec.column_width +"%\">" + $this.list_dat.Records[row].evtriggersec + "</td>\
                            <td class=\"event_list_cell evtriggerele\" style=\"width:"+ $this.column_info.evtriggerele.column_width +"%\">" + $this.list_dat.Records[row].evtriggerele + "</td>\
                            <td class=\"event_list_cell evinvolvedsec\" style=\"width:"+ $this.column_info.evinvolvedsec.column_width +"%\">" + $this.list_dat.Records[row].evinvolvedsec + "</td>\
                            <td class=\"event_list_cell evinvolvedele\" style=\"width:"+ $this.column_info.evinvolvedele.column_width +"%\">" + $this.list_dat.Records[row].evinvolvedele + "</td>\
                            <td class=\"event_list_cell evaveslip\" style=\"width:"+ $this.column_info.evaveslip.column_width +"%\">" + $this.list_dat.Records[row].evaveslip + "</td>\
                            <td class=\"event_list_cell evdetail\" style=\"width:"+ $this.column_info.evdetail.column_width +"%\"></td>\
                        </tr>";
    }
    
    content_str += "</table>";
    
    document.getElementById("event_list").innerHTML = content_str;
    
    $("#event_list tr").mouseenter({$this:$this}, $this.listItemHovered);
    $("#event_list tr").mouseleave({$this:$this}, $this.listItemUnHovered);
    $("#event_list tr").click({$this:$this}, $this.listItemClicked);
    
    $("#event_list tr td.evdetail").mouseenter({$this:$this}, $this.eventDetailHovered);
    $("#event_list tr td.evdetail").mouseleave({$this:$this}, $this.eventDetailUnHovered);
    $("#event_list tr td.evdetail").click({$this:$this}, $this.eventDetailClicked);
},

resizeWindow : function(event) {
    // the FaultMap object
    if ( event == undefined ) {
        $this = this
    } else {
        $this = event.data.$this
    }
    
    the_width = $(window).outerWidth();
    
    content_height = $(window).height() - $("#system").outerHeight(true) - $("#menu_bar").outerHeight(true) - parseFloat($("#uber").css("padding-top").split("px")[0]) - $("#event_list_header").outerHeight(true) - $("#event_list_controls").outerHeight(true);
    
    //position the list controls
    selection_jq = $("#event_list_selection");
    display_items_jq = $("#event_list_display_items");
    
    //console.debug($this._el_controls_dom_jq.css("padding-left"));

    
    //console.debug(content_height);
    
    // This is for when the list is visible
    if ( $this._el_dom_jq.is(":visible") ) {
        
        // Set the visible height of the list
        $this._el_table_dom_jq.height(content_height);
        
        // Position the list controls
        selection_jq.css({"position":"absolute", "left":$this._el_controls_dom_jq.css("padding-left")});
        display_items_jq.css({"position":"absolute", "left":the_width - display_items_jq.outerWidth(true) - parseFloat($this._el_controls_dom_jq.css("padding-right").split("px")[0])});
        
        // Create and position the pager
        available_pager_width = the_width - ($this._el_controls_dom_jq.outerWidth() - $this._el_controls_dom_jq.width()) - selection_jq.outerWidth() - display_items_jq.outerWidth() - 200;

        $this.max_pages = Math.floor((available_pager_width - $this.pagerArrowWidth * 4)/$this.pageWidth);

        $this.buildPager();

        pager_jq = $("#event_list_pager");
        pager_jq.css({"left":the_width/2.0 - pager_jq.outerWidth(true)/2.0});
    }
    
    // This is for when the detail is visible
    if ( $this._ed_dom_jq.is(":visible") ) {
        detail_height = content_height + $("#event_list_header").outerHeight(true) + $("#event_list_controls").outerHeight(true);
        $this._ed_dom_jq.height(detail_height);
        $this._ed_dom_jq.width(the_width);
    }
},

listLoadComplete : function($this) {
    
    //set the height of the list window
    $this.resizeWindow();
    $(window).resize({$this:$this}, $this.resizeWindow);
    
    // for testing only
    //$this.hideEventList($this, 0);
    
},

hideEventList : function($this, event_num) {
    
    $this._el_dom_jq.animate({
        left: "-" + $(window).outerWidth() + "px"
    }, $this.slide_duration, function() {
        $this._el_dom_jq.hide();
        $this._ed_dom_jq.show();
        $this.resizeWindow();
        
        $this._ed_dom_jq.load("event_detail_test_dat.html", {event_num:event_num}, function() {
            $this.detailLoadComplete($this);
        });
    });
    
},

showEventList : function($this) {
    
    system_top = $("#system").offset().top;
    
    //console.debug([$("#system").css("position"),$("#system").css("top"),$("#system").css("width")]);
    $this._el_dom_jq.show();
    
    $this._el_dom_jq.css({"top":0});
    $this._ed_dom_jq.css({"top":-$this._el_dom_jq.height()});
    
    $("#system").css({"position":"absolute","top":system_top, "width":$(window).outerWidth()});
    //$("#system").outerHeight(true)
    //console.debug([$this._el_dom_jq.is(":visible"),$this._ed_dom_jq.is(":visible")])
    
    //$this.resizeWindow();
    
    
    
    $this._ed_dom_jq.animate({
            left: $(window).outerWidth() + "px"
        }, $this.slide_duration);
    
    
    $this._el_dom_jq.animate({
        left: "0px"
    }, $this.slide_duration, function() {
        //$this._el_dom_jq.hide();
        $this._ed_dom_jq.hide();
        $this._ed_dom_jq.css({"top":0, "left":0});
        $this._ed_dom_jq.children().remove();
        $("#system").css({"position":"static", "top":"auto", "width":"100%"});
        $this.resizeWindow();
        
        //$this._ed_dom_jq.load("event_detail_test_dat.html", {event_num:event_num}, function() {
        //    $this.detailLoadComplete($this);
        //});
    });
    
    
},

detailLoadComplete : function($this) {
    
    
    $this._ed_dom_jq.prepend("<div id=\"back_to_event_list\"></div>");
    
    back_button = $("#back_to_event_list");
    
    back_button.mouseenter({$this:$this}, $this.backToListHovered);
    back_button.mouseleave({$this:$this}, $this.backToListUnHovered);
    back_button.click({$this:$this}, $this.backToListClicked);
    
    top_items_height = Math.max($("#event_info").outerHeight() + $("#event_time_series").outerHeight() + back_button.outerHeight(), $("#event_fault_map").outerHeight() + back_button.outerHeight());
    
    $("#event_fault_map").css({"position":"absolute", "top":back_button.outerHeight() + parseFloat(back_button.css("margin-bottom").split("px")[0]), "left":$("#event_time_series").outerWidth() + parseFloat($("#event_fault_map").css("margin-left").split("px")[0])});
    
    $("#event_rupture_map").css({"position":"absolute", "top":top_items_height + parseFloat($("#event_rupture_map").css("margin-top").split("px")[0]) + parseFloat(back_button.css("margin-bottom").split("px")[0])});
    
},

backToListClicked : function(event) {
    // the EventList object
    $this = event.data.$this
    
    $this.showEventList($this);
},

backToListHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    //if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).addClass("hovered");
    //}
    $(event.delegateTarget).css("cursor","pointer");
},

backToListUnHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    //if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).removeClass("hovered");
    //}
    $(event.delegateTarget).css("cursor","auto");
},

eventDetailClicked : function(event) {
    $("html").triggerHandler("click");
    event.stopPropagation();
    
    // the EventList object
    $this = event.data.$this
    
    event_num = parseInt($(event.delegateTarget).parent("tr").children(".column1").text());
    
    $this.hideEventList($this, event_num);
},

eventDetailHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    //if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).addClass("hovered");
    //}
    $(event.delegateTarget).css("cursor","pointer");
},

eventDetailUnHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    //if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).removeClass("hovered");
    //}
    $(event.delegateTarget).css("cursor","auto");
},


listItemHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).addClass("hovered");
    }
    $(event.delegateTarget).css("cursor","pointer");
},

listItemUnHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    if ( !$(event.delegateTarget).hasClass("selected") ) {
        $(event.delegateTarget).removeClass("hovered");
    }
    $(event.delegateTarget).css("cursor","auto");
},

listItemClicked : function(event) {
    // the EventList object
    $this = event.data.$this
    
    event_num = parseInt($(event.delegateTarget).children(".evid").text());
    
    if ( event.metaKey ) {
        if ( $(event.delegateTarget).hasClass("selected") ) {
            if ( $this.server_data.selected[event_num] ) {
                delete $this.server_data.selected[event_num]
            }
            $this.server_data.unselected[event_num] = true;
            $(event.delegateTarget).removeClass("selected");
            //$(event.delegateTarget).find(".event_detail_arrow").attr("src", "images/event_detail_arrow_dark.png");
        } else {
            if ( $this.server_data.unselected[event_num] ) {
                delete $this.server_data.unselected[event_num]
            }
            $this.server_data.selected[event_num] = true;
            $(event.delegateTarget).addClass("selected");
            //$(event.delegateTarget).find(".event_detail_arrow").attr("src", "images/event_detail_arrow_light.png");
        }
    
    } else {
        //console.debug($("tr.selected"))
        $("tr.selected").removeClass("selected").removeClass("hovered");
        $this.server_data.unselected = {all : true};
        $(event.delegateTarget).addClass("selected");
        //$(event.delegateTarget).find(".event_detail_arrow").attr("src", "images/event_detail_arrow_light.png");
        $this.server_data.selected = {};
        $this.server_data.selected[event_num] = true;
    }
    
    console.debug([$this.server_data.selected, $this.server_data.unselected])
    
},

listHeaderHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    if ( !$(event.delegateTarget).hasClass("ascending") || !$(event.delegateTarget).hasClass("descending") ) {
        $(event.delegateTarget).addClass("hovered");
    }
    $(event.delegateTarget).css("cursor","pointer");
},

listHeaderUnHovered : function(event) {
    // the EventList object
    $this = event.data.$this
    
    if ( !$(event.delegateTarget).hasClass("ascending") || !$(event.delegateTarget).hasClass("descending") ) {
        $(event.delegateTarget).removeClass("hovered");
    }
    $(event.delegateTarget).css("cursor","auto");
},

listHeaderClicked : function(event) {
    // the EventList object
    $this = event.data.$this
    
    // remove the current sort
    //console.debug(
    for ( column in $this.column_info ) {
        if ( $this.column_info[column].sort != undefined ) {
            $this.column_info[column].sort = undefined;
            current_sort_column = column;
        }
    }
    
    if ( $(event.delegateTarget).hasClass("ascending") ) {
        sort_direction = "descending";
    } else if ( $(event.delegateTarget).hasClass("descending") ) {
        sort_direction = "ascending";
    } else {
        sort_direction = "ascending";
    }

    $("#event_list_header tr").children("."+current_sort_column).removeClass("ascending descending");
    
    // add the new sort
    if ( $(event.delegateTarget).hasClass("evid") ) {
        selected_column = "evid";
    } else if ( $(event.delegateTarget).hasClass("evyear") ) {
        selected_column = "evyear";
    } else if ( $(event.delegateTarget).hasClass("evmag") ) {
        selected_column = "evmag";
    } else if ( $(event.delegateTarget).hasClass("evtriggersec") ) {
        selected_column = "evtriggersec";
    } else if ( $(event.delegateTarget).hasClass("evtriggerele") ) {
        selected_column = "evtriggerele";
    } else if ( $(event.delegateTarget).hasClass("evinvolvedsec") ) {
        selected_column = "evinvolvedsec";
    } else if ( $(event.delegateTarget).hasClass("evinvolvedele") ) {
        selected_column = "evinvolvedele";
    } else if ( $(event.delegateTarget).hasClass("evaveslip") ) {
        selected_column = "evaveslip";
    }
    
    //console.debug(sort_direction);
    
    $("#event_list_header tr").children("."+selected_column).addClass(sort_direction);
    
    $this.column_info[selected_column].sort = sort_direction;
    
    console.debug("event_list_test_"+$this.number_items_showing+"_dat_" + selected_column + "_" + sort_direction + ".json");
    
    $.getJSON("event_list_test_"+$this.number_items_showing+"_dat_" + selected_column + "_" + sort_direction + ".json", function(data) {
        $this.list_dat = data;
        $this.setPageWidth();
        $this.buildList();
        $this.resizeWindow();
    });
    
}


});