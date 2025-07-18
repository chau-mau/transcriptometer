#!/bin/bash

# transcriptometer.sh
# Author: Dr. Nidhi Sukhija
# Description: Summarizes and visualizes transcriptome/genome assembly statistics from a FASTA file
# Generates key metrics and graphs using gnuplot
# Usage: ./transcriptometer.sh <assembly.fasta>

# Exit immediately on error
set -e

# ========== 0. Input Validation ==========
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <assembly.fasta>"
    exit 1
fi

FASTA=$1

# Check if the input file exists
if [ ! -f "$FASTA" ]; then
    echo "Error: File $FASTA not found!"
    exit 1
fi

# Create output directory for reports
OUTPUT_DIR="./assembly_metrics"
mkdir -p "$OUTPUT_DIR"

echo "===== Assembly Analysis: $FASTA ====="
echo "-------------------------------------"

# ========== 1. Basic Assembly Statistics ==========

echo "1. GENERAL ASSEMBLY STATISTICS:"

# Count number of sequences
TOTAL_SEQS=$(grep -c ">" "$FASTA")

# Total base pairs (excluding headers)
TOTAL_BP=$(grep -v ">" "$FASTA" | tr -d '\n' | wc -c)

# Calculate individual sequence lengths
LENGTHS=$(awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' "$FASTA")

# Min, Max, Average, and Median length
MIN_LEN=$(echo "$LENGTHS" | sort -n | head -1)
MAX_LEN=$(echo "$LENGTHS" | sort -n | tail -1)
AVG_LEN=$(echo "$LENGTHS" | awk '{sum+=$1} END {printf "%.2f", sum/NR}')
MEDIAN_LEN=$(echo "$LENGTHS" | sort -n | awk '{a[NR]=$1} END {print (NR%2==1)?a[int(NR/2)+1]:(a[NR/2]+a[NR/2+1])/2}')

echo "  - Total sequences        : $TOTAL_SEQS"
echo "  - Total length (bp)      : $TOTAL_BP"
echo "  - Min/Avg/Median/Max     : $MIN_LEN / $AVG_LEN / $MEDIAN_LEN / $MAX_LEN"
echo ""

# ========== 2. Nx Statistics ==========

echo "2. Nx STATISTICS:"

# Save lengths in descending order
echo "$LENGTHS" | sort -nr > "$OUTPUT_DIR/seq_lens.txt"

# Compute total base pairs
TOTAL_BP_N50=$(echo "$LENGTHS" | awk '{sum+=$1} END {print sum}')

# Calculate N50 and L50
N50=$(awk -v total="$TOTAL_BP_N50" '{sum+=$1; if (sum >= total*0.5) {print $1; exit}}' "$OUTPUT_DIR/seq_lens.txt")
L50=$(awk -v total="$TOTAL_BP_N50" '{sum+=$1; if (sum >= total*0.5) {print NR; exit}}' "$OUTPUT_DIR/seq_lens.txt")

# Calculate N90 and L90
N90=$(awk -v total="$TOTAL_BP_N50" '{sum+=$1; if (sum >= total*0.9) {print $1; exit}}' "$OUTPUT_DIR/seq_lens.txt")
L90=$(awk -v total="$TOTAL_BP_N50" '{sum+=$1; if (sum >= total*0.9) {print NR; exit}}' "$OUTPUT_DIR/seq_lens.txt")

echo "  - N50                    : $N50 (L50: $L50)"
echo "  - N90                    : $N90 (L90: $L90)"
echo ""

# ========== 3. Base Composition ==========

echo "3. BASE COMPOSITION:"

# GC and AT content as percentages
GC_CONTENT=$(grep -v ">" "$FASTA" | tr -d '\n' | awk '{gc=gsub(/[GCgc]/,""); at=gsub(/[ATat]/,""); printf "%.2f", (gc/(gc+at))*100}')
AT_CONTENT=$(grep -v ">" "$FASTA" | tr -d '\n' | awk '{gc=gsub(/[GCgc]/,""); at=gsub(/[ATat]/,""); printf "%.2f", (at/(gc+at))*100}')

# GC and AT skew: (G-C)/(G+C), (A-T)/(A+T)
GC_SKEW=$(grep -v ">" "$FASTA" | tr -d '\n' | awk '{g=gsub(/[Gg]/,""); c=gsub(/[Cc]/,""); printf "%.3f", (g-c)/(g+c+0.0001)}')
AT_SKEW=$(grep -v ">" "$FASTA" | tr -d '\n' | awk '{a=gsub(/[Aa]/,""); t=gsub(/[Tt]/,""); printf "%.3f", (a-t)/(a+t+0.0001)}')

# Count ambiguous bases (N)
AMBIGUOUS=$(grep -v ">" "$FASTA" | tr -d '\n' | awk '{print gsub(/[Nn]/,"")}')

echo "  - GC Content             : $GC_CONTENT%"
echo "  - AT Content             : $AT_CONTENT%"
echo "  - GC Skew (G-C)/(G+C)    : $GC_SKEW"
echo "  - AT Skew (A-T)/(A+T)    : $AT_SKEW"
echo "  - Ambiguous bases (N)    : $AMBIGUOUS"
echo ""

# ========== 4. Contig Length Distribution ==========

echo "4. CONTIG LENGTH DISTRIBUTION:"

# Bin contigs by length range and print counts
echo "$LENGTHS" | sort -n | awk '
    {if ($1 >= 10000) a["10k+"]++; 
     else if ($1 >= 5000) a["5k-10k"]++; 
     else if ($1 >= 1000) a["1k-5k"]++; 
     else if ($1 >= 500) a["500bp-1k"]++; 
     else a["<500bp"]++} 
    END {for (i in a) print "  - " i ": " a[i] " contigs"}' | sort -k3rn

# ========== 5. Graph Generation ==========

echo ""
echo "5. GENERATING GRAPHS..."
echo "-------------------------------------"

# (A) Length Distribution Histogram
echo "$LENGTHS" | sort -n > "$OUTPUT_DIR/lengths.txt"

gnuplot << EOF
set terminal pngcairo enhanced font "Arial,12"
set output "$OUTPUT_DIR/length_distribution.png"
set title "Sequence Length Distribution (log scale)"
set xlabel "Length (bp)"
set ylabel "Count"
set logscale x
set style fill solid border -1
set grid
binwidth=100
bin(x,width)=width*floor(x/width)
plot "$OUTPUT_DIR/lengths.txt" using (bin(\$1,binwidth)):(1.0) smooth freq with boxes lc rgb "#4CAF50" notitle
EOF

# (B) Nx Curve: Cumulative length over number of contigs
awk '{print $1}' "$OUTPUT_DIR/seq_lens.txt" | sort -nr | awk '{sum+=$1; print NR, sum}' > "$OUTPUT_DIR/nx_data.txt"

gnuplot << EOF
set terminal pngcairo enhanced font "Arial,12"
set output "$OUTPUT_DIR/nx_curve.png"
set title "Nx Curve (Cumulative Assembly)"
set xlabel "Number of Contigs (log scale)"
set ylabel "Cumulative Length (bp)"
set logscale x
set grid
plot "$OUTPUT_DIR/nx_data.txt" using 1:2 with lines lw 2 lc rgb "#2196F3" title "Cumulative Length"
EOF

# (C) GC Content Distribution per contig
awk '/^>/ {if (seq){print gc/len}; seq=1; gc=0; len=0; next} {gc+=gsub(/[GCgc]/,""); len+=length($0)} END {print gc/len}' "$FASTA" > "$OUTPUT_DIR/gc_values.txt"

gnuplot << EOF
set terminal pngcairo enhanced font "Arial,12"
set output "$OUTPUT_DIR/gc_distribution.png"
set title "GC Content Distribution"
set xlabel "GC Content (%)"
set ylabel "Frequency"
set style fill solid border -1
set grid
binwidth=0.02
bin(x,width)=width*floor(x/width)
plot "$OUTPUT_DIR/gc_values.txt" using (bin(\$1,binwidth)*100):(1.0) smooth freq with boxes lc rgb "#FF5722" notitle
EOF

# ========== 6. Summary ==========
echo ""
echo "===== ANALYSIS COMPLETE ====="
echo "Results saved to: $OUTPUT_DIR/"
echo "-------------------------------------"
echo "Graphical Reports:"
echo "  - Length distribution : $OUTPUT_DIR/length_distribution.png"
echo "  - Nx curve            : $OUTPUT_DIR/nx_curve.png"
echo "  - GC distribution     : $OUTPUT_DIR/gc_distribution.png"

# ========== 7. Open Plots (Optional) ==========
# Open plots automatically if GUI tools available
if command -v xdg-open &>/dev/null; then
    xdg-open "$OUTPUT_DIR/length_distribution.png" 2>/dev/null &
    xdg-open "$OUTPUT_DIR/nx_curve.png" 2>/dev/null &
elif command -v open &>/dev/null; then
    open "$OUTPUT_DIR/length_distribution.png" 2>/dev/null &
    open "$OUTPUT_DIR/nx_curve.png" 2>/dev/null &
fi