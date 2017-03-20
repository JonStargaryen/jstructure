function printJSAV(divId, sequences, options)
{
    // Deal with options
    if(options                     == undefined) { options                     = Array();                   }
    if(options.width               == undefined) { options.width               = "400px";                   }
    if(options.height              == undefined) { options.height              = "6pt";                     }
    if(options.nocolor)                          { options.nocolour            = true;                      }
    if(options.toggleNocolor)                    { options.toggleNocolour      = true;                      }
    if(options.fastaLabel          == undefined) { options.fastaLabel          = "Export Selected";         }
    if(options.colorScheme)                      { options.colourScheme        = options.colorScheme;       }
    if(options.colourScheme        == undefined) { options.colourScheme        = "taylor";                  }
    if(options.selectColor)                      { options.selectColour        = true;                      }
    if(options.colorChoices        != undefined) { options.colourChoices       = options.colorChoices;      }
    if(options.deletable)                        { options.selectable          = true;                      }
    if(options.idSubmitAttribute   == undefined) { options.idSubmitAttribute   = "sequence";                }
    if(options.toggleDotifyLabel   == undefined) { options.toggleDotifyLabel   = "Dotify";                  }
    if(options.toggleNocolourLabel == undefined) { options.toggleNocolourLabel = "Do not colour dots";      }
    if(options.toggleNocolorLabel  != undefined) { options.toggleNocolourLabel = options.toggleNocolorLabel;}
    if(options.deleteLabel         == undefined) { options.deleteLabel         = "Delete Selected";         }
    if(options.autoLabels)                       { options.labels              = JSAV_autoLabels(sequences);}

    // Initialize globals if not yet done
    JSAV_init();

    gOptions[divId]         = options;
    gSequences[divId]       = sequences;
    gSequenceLengths[divId] = sequences[0].sequence.length;

    var div = '';
    if($("#" + divId).length == 0) {
        div = $('<div />').appendTo($('body'));
        div.attr('id', divId);
    }
    else{
        div = $("#" + divId);
    }
    var div_sortable = $('<div />').appendTo(div);
    div_sortable.attr('id', divId + '_sortable');
    div_sortable.attr('class', 'JSAVDisplay');

    // 24.09.15
    if(options.scrollY != null)
    {
        div_sortable.css('overflow-y', 'scroll');
        div_sortable.css('white-space', 'nowrap');
        div_sortable.css('height', options.scrollY);
    }

    var html = JSAV_buildSequencesHTML(divId, sequences, options.sortable,
        options.selectable, options.highlight,
        options.dotify, options.nocolour, options.consensus,
        options.labels);
    div_sortable.append(html);


    var div_controls = $('<div />').appendTo(div);
    div_controls.attr('id', divId + '_controls');
    div_controls.attr('class', 'JSAVControls');
}


// ---------------------------------------------------------------------
/**
 Builds the HTML for printing a sequence as a table row. The row
 starts with the identifier and is followed by each amino acid in a
 separate <td> tag with a class to indicate the amino acid type
 (e.g. taylorW for a tryptophan in Willie Taylor scheme).

 @param {string}   id            The identifier
 @param {string}   sequence      A string containing the sequence
 @param {string}   prevSequence  A string containing the previous sequence
 @param {bool}     selectable    Display a selection checkbox
 @param {bool}     dotify        Dotify the sequence
 @param {bool}     nocolour      Don't colour dotified residues
 @param {bool}     isConsensus   This is the consensus sequence
 @param {string}   colourScheme  Name of colour scheme - maps to CSS -
 see JSAV_initColourChoices()
 @param {string}   idSubmit      URL to visit when sequence label clicked
 @param {bool}     idSubmitClean Remove non-alpha characters from sequence
 before submitting it
 @returns {string} text          HTML snippet

 @author
 - 30.05.14 Original  By: ACRM
 - 16.06.14 Added dotify and nocolour - now takes prevSequence parameter
 - 17.06.14 Added isConsensus and colourScheme
 - 18.06.14 Added tooltip
 - 23.09.15 Added idSubmit/idSubmitClean
 - 11.02.16 Added idSubmitAttribute, now takes a sequence object rather than
 the sequence and the id
 */
function JSAV_buildASequenceHTML(sequenceObject, id, sequence, prevSequence, selectable, dotify, nocolour,
                                 isConsensus, colourScheme, idSubmit, idSubmitClean,
                                 idSubmitAttribute)
{
    var seqArray     = sequence.split("");

    var tableLine = "<tr id='" + id + "'>";

    tableLine += "<th class='idCell'>" + id + "</th>";

    var nResidues = seqArray.length;
        for(var i=0; i<nResidues; i++)
        {
            var aa = seqArray[i];
            tableLine += "<td class='" + colourScheme + aa.toUpperCase() + "'>" + aa + "</td>";
        }

    tableLine += "</td></tr>";

    return(tableLine);
}

// ---------------------------------------------------------------------
/**
 Takes an array of sequence objects and builds the HTML to display
 them as a coloured table

 @param   {string}     divId       ID of div in which to print
 @param   {object[]}   sequences   Array of sequence objects
 @param   {bool}       sortable    Should the marker line be displayed
 for sortable displays
 @param   {bool}       selectable  Should check marks be displayed
 @param   {int[]}      highlight   Ranges to be highlighted
 @param   {bool}       dotify      Dotify the sequence alignment
 @param   {bool}       nocolour    Don't colour dotified residues
 @param   {bool}       consensus   Display the consensus sequence
 @param   {array}      labels      Labels to display over sequence
 @returns {string}                 HTML

 @author
 - 30.05.14 Original  By: ACRM
 - 06.06.14 Added call to build the marker row of selected residues
 - 10.06.14 Added sortable and selectable parameters
 - 13.06.14 Added highlight
 - 16.06.14 Added dotify
 - 17.06.14 Added consensus
 - 22.12.15 Added labels
 */
function JSAV_buildSequencesHTML(divId, sequences, sortable, selectable, highlight,
                                 dotify, nocolour, consensus, labels)
{
    var html = "";
    html += "<div class='JSAV'>\n";
    html += "<table border='0'>\n";

    if(labels != undefined)
    {
        html += JSAV_buildLabelsHTML(divId,  gSequenceLengths[divId], selectable, labels);
    }

    // Build the actual sequence entries
    for(var i=0; i<sequences.length; i++)
    {
        var prevSequence = undefined;
        if(i>0) { prevSequence = sequences[i-1].sequence; }
        html += JSAV_buildASequenceHTML(sequences[i], sequences[i].id, sequences[i].sequence, prevSequence,
                selectable, dotify, nocolour, false,
                gOptions[divId].colourScheme,
                gOptions[divId].idSubmit,
                gOptions[divId].idSubmitClean,
                gOptions[divId].idSubmitAttribute) + "\n";
    }

    html += "</table>\n";
    html += "</div>\n";
    return(html);
}


// ---------------------------------------------------------------------
/**
 Initializes global arrays

 @author
 - 09.06.14 Original   By: ACRM
 - 12.06.14 Added more arrays
 - 17.06.14 Added gConsensus array
 */
function JSAV_init()
{
    // Indexed by divId and used to store the values
    try
    {
        if(gSequences == undefined)
        {
            ;
        }
    }
    catch(err)
    {
        gSequences       = Array();
        gOptions         = Array();
        gStartPos        = Array();
        gStopPos         = Array();
        gConsensus       = Array();
        gSequenceLengths = Array();
    }
}




// ---------------------------------------------------------------------
/**
 Create the HTML for a label row in the sequence display

 @param   {string}   divId      - The div we are working in
 @param   {int}      seqLen     - The length of the alignment
 @param   {bool}     selectable - Are there sequences selection boxes
 @param   {string[]} labels     - Array of labels
 @returns {string}              - HTML

 @author
 - 22.12.15 Original   By: ACRM
 */
function JSAV_buildLabelsHTML(divId,  seqLen, selectable, labels)
{
    var html = "";
    if(selectable)
    {
        html += "<tr class='highlightrow'><th></th>";
        html += "<td></td>";
    }
    else
    {
        html += "<tr class='highlightrow'><td></td>";
    }

    for(var i=0; i<labels.length; i++)
    {
        // Make a copy of the label and remove the chain label
        var labelText = labels[i];
        labelText.replace(/^[A-Za-z]/g, '');

        // Find the last character
        var lastChar = labelText.substring(labelText.length-1,labelText.length);

        // Open a table cell with the label as a tooltip
        html += "<td class='tooltip' title='" + labels[i] + "'>";

        // Insert the appropriate character
        if(lastChar == "0")                   // 0 - do a '|'
        {
            html += "|";
        }
        else if (lastChar.match(/[A-Za-z]/))  // Insert code - show the code
        {
            html += lastChar;
        }
        else                                  // Otherwise do a '.'
        {
            html += ".";
        }

        // And finish the table cell
        html += "</td>";
    }

    // Finish the table row
    html += "</tr>\n";

    return html;
}


// ---------------------------------------------------------------------
/**
 Create an array of labels

 @param   {object[]} sequences  - Array of sequence objects
 @returns {string[]}            - Array of labels

 @author
 - 22.12.15 Original   By: ACRM
 */
function JSAV_autoLabels(sequences)
{
    var seqLen = sequences[0].sequence.length;
    var labels = Array();
    for(var i=1; i<=seqLen; i++)
    {
        labels.push(i.toString());
    }
    return labels;
}