proc _list_arithmetic {l op} { 
    for {set i 0} {$i < [llength $op]} {incr i} {
        if { [string is integer [lindex $op $i]]} {
            set op [lreplace $op $i $i [lindex $l [lindex $op $i]]] 
        } 
    }
    return [expr {*}$op]
}

proc load_csv_data {fileName col {targetField beta} {bIgnoreAtom 1} {mol top}} {
    # The target data is a pandas DataFRame output where the first four columns are reserved for identifiation purposes.
    # Nominally these are: key, segname, resid, name
    # Note: This is zero-indexed so the first data column is 0.
    # Secondary usage is to sum and subtract columns. Will parse columns in the format
    # 2-1 => column 2 minus column 1.
    # Another secondary usage is to put in a pure alpha as a $col argument.
    # This performs an lsearch on the headers to obtain the fieldID.
    # Note that it should not be amongst one of the first 4 columns!

    set fp [open $fileName r]
    set dat [split [read $fp] "\n"]
    close $fp
    set bFirst 1
    foreach line $dat {
        if { [llength $line] == 0} {continue}
        if { [lindex $line 0] == "#" } {continue}
        set line [split $line ","]
        if { $bFirst } {
            set bFirst 0
            set selType1 [lindex $line 1]
            set selType2 [lindex $line 2]
            set selType3 [lindex $line 3]
            if {[string is alpha $col]} {
                set col [expr [lsearch $line $col]-4 ]
                if { [$col < 0] } {
                    puts "= = = ERROR: alphabetical column argument $col not found in $line or is one of the reserved identification columns!"
                    return
                } 
            }
            # = = = Interpret string inputs
            continue
        }
        # = = = First column as ID of residue, atom, etc.
        if { $bIgnoreAtom } {
            set selText "$selType1 [lindex $line 1] and $selType2 [lindex $line 2]"
        } else {
            set selText "$selType1 [lindex $line 1] and $selType2 [lindex $line 2] and $selType3 [lindex $line 3]"
        }
        set tempSelection [atomselect $mol $selText]
        if {[$tempSelection num] == 0} {
            puts "= = = WARNING: no atoms selected by $selText! "
            return
        }
        set line [lreplace $line 0 3]
        if { [llength $col] > 1 } {
            $tempSelection set $targetField [_list_arithmetic $line $col]
        } else {
            $tempSelection set $targetField [lindex $line $col]
        }
        $tempSelection delete
    }
    puts "= = = CSV data loaded from $fileName"
}
