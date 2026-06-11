from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, PageBreak, HRFlowable)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY, TA_RIGHT

OUTFILE = "/Users/boris/Desktop/RSShapley/Combination_Monte/Final/codess/JMP/RSM_summary.pdf"

# ── Palette ───────────────────────────────────────────────────────────────────
BLUE  = colors.HexColor("#1a4f8a")
LBLUE = colors.HexColor("#e8eef6")
MGREY = colors.HexColor("#dce3ec")
DGREY = colors.HexColor("#6c757d")
GREEN = colors.HexColor("#2e7d32")
RED   = colors.HexColor("#c62828")
AMBER = colors.HexColor("#e65100")
WHITE = colors.white; BLACK = colors.black

doc = SimpleDocTemplate(OUTFILE, pagesize=letter,
    leftMargin=0.8*inch, rightMargin=0.8*inch,
    topMargin=0.85*inch, bottomMargin=0.85*inch)

styles = getSampleStyleSheet()
def S(name, parent="Normal", **kw):
    return ParagraphStyle(name, parent=styles[parent], **kw)

Title   = S("T", "Title",   fontSize=18, textColor=BLUE, spaceAfter=3, alignment=TA_CENTER)
Sub     = S("Su","Normal",  fontSize=9.5, textColor=DGREY, spaceAfter=12, alignment=TA_CENTER)
H1      = S("H1","Heading1",fontSize=13, textColor=BLUE, spaceBefore=12, spaceAfter=3)
H2      = S("H2","Heading2",fontSize=10.5, textColor=BLUE, spaceBefore=7, spaceAfter=2)
Body    = S("Bo","Normal",  fontSize=8.8, leading=13, spaceAfter=5, alignment=TA_JUSTIFY)
Caption = S("Ca","Normal",  fontSize=7.8, textColor=DGREY, spaceAfter=5, alignment=TA_CENTER)

# Cell styles for tables
def cs(size=8, bold=False, color=BLACK, align=TA_CENTER):
    return ParagraphStyle("c", fontSize=size, leading=size+2.5,
        fontName="Helvetica-Bold" if bold else "Helvetica",
        textColor=color, alignment=align, wordWrap='CJK')

HDR  = cs(8, bold=True, color=WHITE)
CELL = cs(8)
CELLL= cs(8, align=TA_LEFT)
SMLL = cs(7.5, align=TA_LEFT)
SMLC = cs(7.5, align=TA_CENTER)

def p(text, style=CELL): return Paragraph(text, style)
def pl(text, style=CELLL): return Paragraph(text, style)

def hline(thick=0.5, col=MGREY):
    return HRFlowable(width="100%", thickness=thick, color=col, spaceAfter=5, spaceBefore=1)

BASE_TS = [
    ("BACKGROUND",   (0,0),(-1,0), BLUE),
    ("TEXTCOLOR",    (0,0),(-1,0), WHITE),
    ("FONTNAME",     (0,0),(-1,-1),"Helvetica"),
    ("FONTSIZE",     (0,0),(-1,-1), 8),
    ("ALIGN",        (0,0),(-1,-1),"CENTER"),
    ("VALIGN",       (0,0),(-1,-1),"MIDDLE"),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[WHITE, LBLUE]),
    ("GRID",         (0,0),(-1,-1), 0.35, MGREY),
    ("TOPPADDING",   (0,0),(-1,-1), 3),
    ("BOTTOMPADDING",(0,0),(-1,-1), 3),
    ("LEFTPADDING",  (0,0),(-1,-1), 4),
    ("RIGHTPADDING", (0,0),(-1,-1), 4),
]

def tbl(data, widths, extra=None):
    t = Table(data, colWidths=widths, repeatRows=1)
    ts = TableStyle(BASE_TS)
    if extra:
        for cmd in extra: ts.add(*cmd)
    t.setStyle(ts)
    return t

story = []

# ═══ COVER ════════════════════════════════════════════════════════════════════
story += [Spacer(1, 0.4*inch),
    Paragraph("Random Subset Methods for Forecast Combination", Title),
    Paragraph("Theory &nbsp;·&nbsp; Monte Carlo Simulations &nbsp;·&nbsp; Empirical Results — Alignment Report", Sub),
    Paragraph("4 June 2026", Sub),
    hline(1.5, BLUE), Spacer(1,4),
    Paragraph(
        "This document summarises the current state of the RSM project across three components: "
        "(i) theoretical results derived in the proposal and addendum, "
        "(ii) exogenous and endogenous Monte Carlo simulations, and "
        "(iii) the GDP empirical application using FRED-MD/FRED-QD real-time vintages. "
        "Each theoretical prediction is evaluated against its simulation and empirical counterpart.",
        Body)]

# ═══ 1. THEORY ════════════════════════════════════════════════════════════════
story += [Paragraph("1. Theory", H1), hline(1.5, BLUE)]
story += [Paragraph("1.1 Model", H2),
    Paragraph(
        "The common-loading model specifies <b>&#931;<sub>e</sub> = q&middot;&#953;&#953;' + D</b>, "
        "D = diag(&#963;<sub>1</sub><super>2</super>,...,&#963;<sub>m</sub><super>2</super>), "
        "precision weights &#964;<sub>i</sub> = &#963;<sub>i</sub><super>&minus;2</super>. "
        "Two population risk objects are central:", Body)]

story.append(tbl([
    [p("Object",HDR), p("Formula",HDR), p("Interpretation",HDR)],
    [p("Q<sub>k</sub><super>sub</super>"),
     p("q + E<sub>S</sub>[1/A<sub>S</sub>]"),
     pl("Average-subset risk (B subsets averaged)")],
    [p("Q<sub>k</sub><super>bag,&#8734;</super>"),
     p("v&#772;<sub>k</sub>' &#931;<sub>e</sub> v&#772;<sub>k</sub>"),
     pl("Infinite-bagging risk (G&#8594;&#8734;)")],
    [p("Q<sub>k</sub><super>bag,G</super>"),
     p("Q<super>bag,&#8734;</super> + (1/G)(Q<super>sub</super> &#8722; Q<super>bag,&#8734;</super>)"),
     pl("Finite-bagging bridge — <b>exact identity</b>")],
], [1.3*inch, 2.4*inch, 3.0*inch]))

story += [Spacer(1,6), Paragraph("1.2 Key Theoretical Results", H2)]
story.append(tbl([
    [p("Result",HDR), p("Formula",HDR), p("Assumption",HDR)],
    [pl("Small-het excess risk"),
     p("R<sub>k</sub> &#8776; C<sub>&#948;</sub> z<sub>k</sub><super>2</super>"),
     pl("&#964;<sub>i</sub> &#8776; &#964;&#772;(1+&#948;<sub>i</sub>), &#948;<sub>i</sub> small")],
    [pl("Cube-root optimal k (bagged)"),
     p("k* &#8776; (2C<sub>&#948;</sub>T/A)<super>1/3</super>"),
     pl("Small-het + T&#8811;k + z<sub>k</sub>&#8776;1/k")],
    [pl("Square-root optimal k (subset)"),
     p("k* &#8776; (A<sub>0</sub>T/B)<super>1/2</super>"),
     pl("Leading-order average-subset loss")],
    [pl("Finite-bagging bridge"),
     p("Q<super>bag,G</super> = Q<super>bag,&#8734;</super> + G<super>&#8722;1</super>(Q<super>sub</super>&#8722;Q<super>bag,&#8734;</super>)"),
     pl("Exact algebraic identity, no approximation")],
    [pl("RSM dominance (&#954;)"),
     p("Q<super>bag,&#8734;</super> &lt; L<super>full OLS</super> for k &#8805; &#954;"),
     pl("&#954; = 2 in all simulated DGPs")],
], [1.65*inch, 2.25*inch, 2.85*inch]))

story += [Spacer(1,6), Paragraph("1.3 Power-Law Generalisation (Addendum, Section 13)", H2),
    Paragraph(
        "The small-het formula underestimates the true excess risk. "
        "The addendum proposes <b>R<sub>k</sub> &#8776; C&#183;z<sub>k</sub><super>&#945;</super></b> "
        "and establishes two limit cases:", Body)]

story.append(tbl([
    [p("Regime",HDR), p("True form of R<sub>k</sub>",HDR),
     p("Slope &#945;&#770; (OLS)",HDR), p("k* scaling",HDR)],
    [pl("Homoskedastic"), p("0"), p("—"), p("Any k optimal")],
    [pl("Small heterogeneity"),
     p("C<sub>&#948;</sub> z<sub>k</sub><super>2</super>"),
     p("2 (exact)"), p("T<super>1/3</super>")],
    [pl("Intermediate"), p("Mixture"),
     p("&#945; &#8712; (1, 2)"), p("T<super>1/(&#945;+1)</super>")],
    [pl("Large het (dominant forecaster)"),
     p("C'(kz<sub>k</sub>)<super>2</super>"),
     p("&#8776; 1 (grid OLS)"), p("T<super>1/2</super>")],
], [1.75*inch, 1.85*inch, 1.5*inch, 1.65*inch]))

story += [Spacer(1,4),
    Paragraph(
        "Generalised optimal k* = (C&#945;T/A)<super>1/(&#945;+1)</super>. "
        "A general proof that &#945; &#8712; [1, 2] for arbitrary {&#964;<sub>i</sub>} "
        "<b>remains open</b> (noted explicitly in the addendum).", Body)]

# ═══ 2. EXOGENOUS MC ══════════════════════════════════════════════════════════
story += [PageBreak(), Paragraph("2. Exogenous Monte Carlo", H1), hline(1.5, BLUE),
    Paragraph(
        "Design: m=50, T=240, &#931;<sub>e</sub> <b>known</b> (exogenous weights). "
        "N_REP=500, G=800 bags, N_SUBSETS=2500. "
        "Four DGPs: B &#8712; {3,8} &#215; H &#8712; {1.05, 6} "
        "(H = heterogeneity ratio max/min &#963;<sub>i</sub><super>2</super>). "
        "Status: <b>baseline complete; scaling_T_m complete; joint_B_H_grid running</b>.", Body),
    Paragraph("2.1 Baseline Results (m=50, T=240)", H2)]

story.append(tbl([
    [p("DGP",HDR), p("k*(formula)",HDR), p("k*(grid)",HDR),
     p("k*(empirical)",HDR), p("Gain vs Equal",HDR),
     p("Gap vs Oracle",HDR), p("&#945;&#770;",HDR), p("R<super>2</super>",HDR)],
    [pl("B=3, H=1.05"), p("1"), p("1"),
     p("3",cs(8,color=AMBER)), p("0.04%"), p("0.19%"), p("1.05"), p("0.982")],
    [pl("B=8, H=1.05"), p("1"), p("1"),
     p("3",cs(8,color=AMBER)), p("0.09%"), p("0.45%"), p("1.05"), p("0.976")],
    [pl("B=3, H=6"), p("3",cs(8,bold=True,color=GREEN)), p("3",cs(8,bold=True,color=GREEN)),
     p("6",cs(8,bold=True,color=GREEN)),
     p("10.4%",cs(8,bold=True,color=GREEN)), p("0.92%"), p("1.66"), p("0.986")],
    [pl("B=8, H=6"), p("5",cs(8,bold=True,color=GREEN)), p("5",cs(8,bold=True,color=GREEN)),
     p("9",cs(8,bold=True,color=GREEN)),
     p("42.1%",cs(8,bold=True,color=GREEN)), p("3.23%"), p("1.68"), p("0.983")],
], [1.15*inch, 0.72*inch, 0.65*inch, 0.82*inch, 0.82*inch, 0.82*inch, 0.55*inch, 0.52*inch]))

story += [Spacer(1,3),
    Paragraph("Cube-root formula matches population grid k* exactly in all four DGPs. "
        "Empirical k* &#8776; 2&#215; population k* (&#931; estimation noise penalty). "
        "Gains only meaningful when H is large.", Caption),
    Paragraph("2.2 Small-Het Approximation Quality", H2),
    Paragraph("The ratio R<sub>mc</sub>/R<sub>small-het</sub> grows with k and blows up near k=m. "
        "Counterintuitively, approximation error is <b>worst for H=1.05</b> "
        "(near-homogeneous): the true curve follows C'(kz<sub>k</sub>)<super>2</super> "
        "while the formula gives C<sub>&#948;</sub>z<sub>k</sub><super>2</super>, "
        "diverging at rate k<super>2</super>.", Body)]

story.append(tbl([
    [p("DGP",HDR), p("Ratio at k=5",HDR), p("Ratio at k=20",HDR), p("Mean rel. error",HDR)],
    [pl("B=3, H=1.05"), p("3.9&#215;",cs(8,bold=True,color=RED)),
     p("21&#215;",cs(8,bold=True,color=RED)), p("91%",cs(8,bold=True,color=RED))],
    [pl("B=8, H=1.05"), p("3.3&#215;",cs(8,bold=True,color=RED)),
     p("16&#215;",cs(8,bold=True,color=RED)), p("90%",cs(8,bold=True,color=RED))],
    [pl("B=3, H=6"),    p("2.6&#215;",cs(8,color=AMBER)), p("3.8&#215;",cs(8,color=AMBER)), p("72%",cs(8,color=AMBER))],
    [pl("B=8, H=6"),    p("2.3&#215;",cs(8,color=AMBER)), p("9.3&#215;",cs(8,color=AMBER)), p("78%",cs(8,color=AMBER))],
], [1.5*inch, 1.4*inch, 1.4*inch, 1.45*inch]))

story += [Spacer(1,6), Paragraph("2.3 Scaling Experiment (T &#8712; {80, 240, 640, 1000}, H=1.05)", H2),
    Paragraph("The H=1.05 DGP has near-zero excess risk relative to Q<sub>opt</sub>&#8776;1.06, "
        "so the population k*(grid)=1 for all T — the cube-root rule is not operational. "
        "The H=6 scaling run (500 reps, paper profile) is <b>currently running</b> "
        "and will provide the definitive T-scaling check.", Body)]

story.append(tbl([
    [p("T",HDR), p("k*(pop. grid)",HDR), p("k*(cube-root, cont.)",HDR), p("k*(empirical)",HDR)],
    [p("80"), p("1"), p("0.57"), p("2")],
    [p("240"), p("1"), p("0.82"), p("2")],
    [p("640"), p("1"), p("1.13"), p("4")],
    [p("1000"), p("1"), p("1.31"), p("3")],
], [1.0*inch, 1.3*inch, 1.7*inch, 1.3*inch]))
story.append(Paragraph(
    "Empirical log-log slope: 0.238 (p=0.21, not significant; k* pinned near floor). "
    "H=6 run pending — expected to show clear cube-root exponent.", Caption))

# ═══ 3. ENDOGENOUS MC ═════════════════════════════════════════════════════════
story += [PageBreak(), Paragraph("3. Endogenous Monte Carlo", H1), hline(1.5, BLUE),
    Paragraph("Design: m=50, T=200, &#931;<sub>e</sub> <b>estimated from data</b> (endogenous weights). "
        "Two error structures (Toeplitz &#961;=0.8, Factor r=2) &#215; "
        "two &#946; designs (B2 dense, B4 sparse). N_REP=200, K_RS=200, OOS_LEN=30. "
        "Status: <b>complete</b>.", Body),
    Paragraph("3.1 Four-Scenario Summary", H2)]

story.append(tbl([
    [p("Scenario",HDR), p("k*(RMSFE)",HDR), p("Rel. RMSFE",HDR),
     p("Gain",HDR), p("k*(L<sub>emp</sub>)",HDR), p("&#954;",HDR)],
    [pl("Toeplitz &#961;=0.8, B2 (dense)"),  p("11"), p("0.965"), p("&#8722;3.5%"), p("32"), p("2")],
    [pl("Toeplitz &#961;=0.8, B4 (sparse)"), p("11"), p("0.902",cs(8,bold=True,color=GREEN)),
     p("&#8722;9.8%",cs(8,bold=True,color=GREEN)), p("26"), p("2")],
    [pl("Factor r=2, B2 (dense)"),  p("14"), p("0.950"), p("&#8722;5.0%"), p("29"), p("2")],
    [pl("Factor r=2, B4 (sparse)"), p("14"), p("0.939",cs(8,bold=True,color=GREEN)),
     p("&#8722;6.1%",cs(8,bold=True,color=GREEN)), p("29"), p("2")],
], [2.1*inch, 0.75*inch, 0.85*inch, 0.75*inch, 0.85*inch, 0.5*inch]))
story.append(Paragraph(
    "&#954;=2 universally: bagged RSM beats full estimated OLS for any k&#8805;2. "
    "Sparsity (B4) amplifies gains. Factor structure shifts k* 11&#8594;14.", Caption))

story += [Spacer(1,4), Paragraph("3.2 T-Scaling of Empirical k* (B4, 120 reps, T &#8712; {100,...,800})", H2)]
story.append(tbl([
    [p("&#931; type",HDR), p("T=100",HDR), p("T=200",HDR), p("T=300",HDR),
     p("T=500",HDR), p("T=800",HDR), p("Log-log slope",HDR),
     p("Slope=1/3?",HDR), p("Slope=1/2?",HDR)],
    [pl("Toeplitz"), p("10"), p("14"), p("14"), p("18"), p("22"),
     p("0.394**",cs(8,bold=True,color=GREEN)), p("Cannot reject"), p("Rejected (t=&#8722;2.1)")],
    [pl("Factor"), p("14"), p("14"), p("18"), p("26"), p("26"),
     p("0.367**",cs(8,bold=True,color=GREEN)), p("Cannot reject"), p("Rejected (t=&#8722;1.9)")],
], [0.75*inch,0.42*inch,0.42*inch,0.42*inch,0.42*inch,0.42*inch,0.9*inch,0.9*inch,1.05*inch]))
story.append(Paragraph(
    "Log-log slope (0.37&#8211;0.39) is statistically indistinguishable from 1/3 but significantly below 1/2. "
    "Appearance of k*&#8776;&#8730;T at T=100&#8211;200 is a coincidence: "
    "c&#183;T<super>1/3</super> = &#8730;T when T = c<super>6</super> &#8776; 113 (c&#8776;2.2). "
    "At T=800: cube-root gives &#8776;20 (matches data=22); &#8730;T gives 28.3 (does not).", Body))

story += [Spacer(1,4), Paragraph("3.3 Population vs Empirical k* Gap", H2),
    Paragraph("The endogenous simulation reveals a systematic three-way split "
        "absent in the exogenous case:", Body)]
story.append(tbl([
    [p("Criterion",HDR), p("Loss surface",HDR), p("k* (Toeplitz B2)",HDR)],
    [pl("L<sub>pop</sub> (known &#931;<sub>e</sub>)"),
     p("Range 0.009 — nearly flat",cs(8,color=AMBER)), p("2",cs(8,color=AMBER))],
    [pl("RMSFE (OOS, empirical)"), p("—"), p("11",cs(8,bold=True))],
    [pl("L<sub>emp</sub> (estimated &#931;&#770;<sub>e</sub>)"),
     p("Range 0.39 — steep",cs(8,color=RED)), p("32",cs(8,bold=True,color=RED))],
], [2.2*inch, 2.4*inch, 1.6*inch]))
story.append(Paragraph(
    "The inflation formula (T&#8722;1)/(T&#8722;k) accounts for estimating k combination weights given "
    "fixed &#931;. It does <b>not</b> account for estimating &#931; from m=50 series — "
    "this additional noise pushes empirical k* from 2 to 11.", Body))

# ═══ 4. GDP EMPIRICAL ═════════════════════════════════════════════════════════
story += [PageBreak(), Paragraph("4. GDP Empirical Application (FRED-MD/FRED-QD)", H1),
    hline(1.5, BLUE),
    Paragraph("102 real-time vintages (1999Q3&#8211;2025Q1, vintage 18 excluded). "
        "Predictor pool: m&#8776;126 monthly FRED-MD series, each transformed to a "
        "predicted quarterly growth rate by a rolling univariate AR (n<sub>start</sub>=40 burn-in). "
        "k-grid: 2&#8211;55. RSM draws: B<sub>main</sub>=1000, B<sub>allk</sub>=200. "
        "Status: <b>complete</b> (minor legend bug in combined figure; all data saved).", Body),
    Paragraph("4.1 Sample Design: Why T&#8776;81 and m&#8804;55", H2)]

story.append(tbl([
    [p("Parameter",HDR), p("Value",HDR), p("Source / explanation",HDR)],
    [pl("Training obs T"), pl("79&#8211;180 (median 81)"),
     pl("FRED-QD starts &#8776;1970Q1; after n<sub>start</sub>=40 burn-in: "
        "T &#8776; T<sub>QD</sub> &#8722; 40")],
    [pl("Rule for k<sub>main</sub>"), pl("round(&#8730;T) = 9 (median)"),
     pl("Hard-coded heuristic; cube-root formula <b>not used</b>")],
    [pl("Predictor pool m"), pl("&#8776;126 FRED-MD monthly series"),
     pl("Monthly-to-quarterly aggregation; each is a univariate AR forecast")],
    [pl("k search range"), pl("2 to 55"),
     pl("Hard-coded cap; at k=55, T=81: k/T&#8776;0.68 — deeply crowded OLS, "
        "inflation factor (T&#8722;1)/(T&#8722;k)&#8776;3.1")],
    [pl("C<sub>&#948;</sub> estimated?"), pl("No"),
     pl("m&gt;T makes &#931;<sub>e</sub> estimation infeasible; "
        "C<sub>&#948;</sub>, &#945;, A all unobservable")],
], [1.1*inch, 1.45*inch, 4.2*inch]))

story += [Spacer(1,6), Paragraph("4.2 Forecast Accuracy Tables", H2),
    Paragraph("<b>Table 1 — Full Sample</b>", Body)]
story.append(tbl([
    [p("Model",HDR), p("MAFE",HDR), p("Sig.",HDR), p("RMSFE",HDR), p("Sig.",HDR),
     p("Target MAFE",HDR), p("Target RMSFE",HDR)],
    [pl("Mean"),          p("0.54"), p(""),   p("1.10"), p(""),   p("0.54"), p("1.10")],
    [pl("RSM"),           p("0.45"), p("**"), p("0.76"), p("*"),  p("0.45"), p("0.74")],
    [pl("Lasso"),         p("0.51"), p(""),   p("0.73"), p("*"),  p("0.51"), p("0.76")],
    [pl("Ridge"),         p("0.49"), p("*"),  p("0.98"), p(""),   p("0.50"), p("1.00")],
    [pl("Random Forest"), p("0.55"), p(""),   p("1.24"), p(""),   p("0.55"), p("1.24")],
], [1.15*inch,0.65*inch,0.45*inch,0.65*inch,0.45*inch,0.9*inch,0.92*inch]))
story.append(Spacer(1,5))
story.append(Paragraph("<b>Table 2 — COVID-19 Excluded</b>", Body))
story.append(tbl([
    [p("Model",HDR), p("MAFE",HDR), p("Sig.",HDR), p("RMSFE",HDR), p("Sig.",HDR),
     p("Target MAFE",HDR), p("Target RMSFE",HDR)],
    [pl("Mean"),          p("0.42"), p(""),   p("0.60"), p(""),   p("0.42"), p("0.60")],
    [pl("RSM"),           p("0.38"), p("**"), p("0.53"), p("**"), p("0.39"), p("0.53")],
    [pl("Lasso"),         p("0.45"), p(""),   p("0.58"), p(""),   p("0.45"), p("0.60")],
    [pl("Ridge"),         p("0.39"), p("**"), p("0.54"), p("**"), p("0.39"), p("0.55")],
    [pl("Random Forest"), p("0.41"), p(""),   p("0.57"), p("*"),  p("0.41"), p("0.58")],
], [1.15*inch,0.65*inch,0.45*inch,0.65*inch,0.45*inch,0.9*inch,0.92*inch]))
story.append(Paragraph("RSM matches paper targets within &#177;0.02 RMSFE.", Caption))

story += [Spacer(1,4),
    Paragraph("4.3 Optimal k vs Theory Benchmarks (T&#8776;81)", H2)]
story.append(tbl([
    [p("Benchmark",HDR), p("Full sample",HDR), p("No-COVID",HDR)],
    [pl("k*(RMSFE, actual)"),
     p("19",cs(8,bold=True,color=RED)), p("21",cs(8,bold=True,color=RED))],
    [pl("&#8730;T  (T=81)"),                             p("9.0"), p("9.0")],
    [pl("T<super>1/3</super>  (T=81)"),                  p("4.3"), p("4.3")],
    [pl("2.2 &#215; T<super>1/3</super>  (endo. Toeplitz constant)"), p("9.5"), p("9.5")],
    [pl("2.8 &#215; T<super>1/3</super>  (endo. Factor constant)"),   p("12.1"), p("12.1")],
    [pl("Power-law slope &#945;&#770;  (k=2..30)"),
     p("0.10",cs(8,color=AMBER)), p("0.04",cs(8,color=AMBER))],
], [3.3*inch, 1.25*inch, 1.25*inch]))
story.append(Paragraph(
    "Empirical k*=19&#8211;21 is &#8776;2&#215; larger than endogenous simulation predicts. "
    "The RMSFE curve is nearly flat (&#945;&#770;&#8776;0.04&#8211;0.10, far below theoretical [1,2]); "
    "102 OOS predictions give high sampling variance. The &#8730;T heuristic (k=9) "
    "costs only 0.09&#8211;0.09 RMSFE relative to the empirical optimum.", Body))

# ═══ 5. ALIGNMENT ═════════════════════════════════════════════════════════════
story += [PageBreak(), Paragraph("5. Theory&#8211;Simulation&#8211;Empirical Alignment", H1),
    hline(1.5, BLUE),
    Paragraph("Each row is a theoretical prediction; columns show whether it is "
        "confirmed (&#10003; green), partial (&#8764; amber), not confirmed (&#10007; red), "
        "or not applicable (&#8212;) in each component.", Body)]

CK = cs(7.5, bold=True, color=GREEN)
PA = cs(7.5, bold=True, color=AMBER)
NO = cs(7.5, bold=True, color=RED)
NA = cs(7.5, color=DGREY)

def ck(t): return p(t, CK)
def pa(t): return p(t, PA)
def no(t): return p(t, NO)
def na(t="&#8212; N/A"): return p(t, NA)

rows = [
    [p("Theoretical Prediction",HDR), p("Exogenous MC",HDR),
     p("Endogenous MC",HDR), p("GDP Empirical",HDR)],

    [pl("Finite-bagging bridge\nQ<super>bag,G</super> = Q<super>bag,&#8734;</super> + "
        "(1/G)(Q<super>sub</super>&#8722;Q<super>bag,&#8734;</super>)",SMLL),
     ck("&#10003; Confirmed\nmax|err| &lt; 9e-15"), na(), na()],

    [pl("Cube-root formula matches\npopulation grid k*",SMLL),
     ck("&#10003; Exact\n(all 4 DGPs)"),
     pa("&#8764; Partial\nPop. k*=2; emp. k*=11\n(&#931; estimation noise)"),
     no("&#10007; Not applicable\nC not estimable (m&gt;T)")],

    [pl("Square-root subset formula\nmatches grid k*",SMLL),
     ck("&#10003; Exact\n(all 4 DGPs)"), na(), na()],

    [pl("Empirical k* scales as T<super>1/3</super>",SMLL),
     pa("&#8764; Partial\nH=1.05: k* at floor\nH=6 run pending"),
     ck("&#10003; Confirmed\nSlope 0.37&#8211;0.39 (p&lt;0.01)"),
     pa("&#8764; Partial\nk* grows; curve too flat\nto measure slope")],

    [pl("Power-law &#945;&#770; &#8712; [1,2], R&#178;&gt;0.97",SMLL),
     ck("&#10003; Confirmed\n&#945;&#770; &#8712; {1.05, 1.68}"), na(),
     no("&#10007; &#945;&#770; &#8776; 0.04&#8211;0.10\nsampling noise dominates")],

    [pl("&#945;&#770; determined by H,\nnot by B (signal strength)",SMLL),
     ck("&#10003; Confirmed\nB=3 and B=8 identical"), na(), na()],

    [pl("RSM dominates full\nestimated OLS (&#954;=2)",SMLL),
     ck("&#10003; ~25% lower loss"),
     ck("&#10003; &#954;=2 universally"),
     ck("&#10003; OLS infeasible;\nRSM well-defined")],

    [pl("RSM beats equal weights\nwhen H is large",SMLL),
     ck("&#10003; 10&#8211;42% gains"),
     ck("&#10003; 3.5&#8211;9.8% RMSFE"),
     ck("&#10003; 5&#8211;12% RMSFE gain")],

    [pl("Small-het approximation\naccurate quantitatively",SMLL),
     no("&#10007; Off by 2&#8211;21&#215;\n72&#8211;91% rel. error"), na(), na()],

    [pl("k* &#8776; &#8730;T as rule of thumb",SMLL),
     pa("&#8764; Coincidence at T&#8776;113"),
     pa("&#8764; Valid at T=100&#8211;200\ntrue scaling is T<super>1/3</super>"),
     pa("&#8764; Used in script;\ngrid k* is ~2&#215; larger")],
]

align_tbl = Table(rows, colWidths=[2.1*inch, 1.5*inch, 1.55*inch, 1.6*inch], repeatRows=1)
align_tbl.setStyle(TableStyle([
    ("BACKGROUND",   (0,0),(-1,0), BLUE),
    ("FONTNAME",     (0,0),(-1,-1),"Helvetica"),
    ("FONTSIZE",     (0,0),(-1,-1), 7.5),
    ("ALIGN",        (0,0),(-1,-1),"CENTER"),
    ("VALIGN",       (0,0),(-1,-1),"TOP"),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[WHITE,LBLUE]),
    ("GRID",         (0,0),(-1,-1), 0.35, MGREY),
    ("TOPPADDING",   (0,0),(-1,-1), 3),
    ("BOTTOMPADDING",(0,0),(-1,-1), 3),
    ("LEFTPADDING",  (0,0),(-1,-1), 4),
    ("RIGHTPADDING", (0,0),(-1,-1), 4),
    ("ALIGN",        (0,0),(0,-1),"LEFT"),
]))
story.append(align_tbl)

# ═══ 6. OPEN ITEMS ════════════════════════════════════════════════════════════
story += [PageBreak(), Paragraph("6. Open Items and Next Steps", H1), hline(1.5, BLUE)]

def pri(label, color):
    return p(f"<b>{label}</b>", cs(8, bold=True, color=color))

open_rows = [
    [p("Priority",HDR), p("Item",HDR), p("Status",HDR)],
    [pri("High",RED),
     pl("H=6 scaling run (500 reps, paper profile): verify cube-root "
        "exponent &#8776;1/3 for exogenous case with meaningful heterogeneity",SMLL),
     pl("Running — PID active, &#8764;4&#8211;5 h remaining",SMLL)],
    [pri("High",RED),
     pl("joint_B_H_grid: full &#945;&#770;(H) surface, B&#215;H interaction, "
        "cube-root scaling at H=6",SMLL),
     pl("Running — begins after scaling_T_m finishes",SMLL)],
    [pri("High",RED),
     pl("Prove &#945; &#8712; [1,2] for arbitrary {&#964;<sub>i</sub>}: "
        "upper bound (Jensen argument) plausible; lower bound needs "
        "monotonicity in H or per-forecaster decomposition",SMLL),
     pl("Open theoretical problem",SMLL)],
    [pri("High",RED),
     pl("Explain robustness of cube-root formula: why does k*(formula)=k*(grid) "
        "even when small-het approx fails by 3&#8211;21&#215;?",SMLL),
     pl("Empirical finding — no proof yet",SMLL)],
    [pri("Medium",AMBER),
     pl("Endogenous G_bag fix: re-run with G_bag=K_RS to populate "
        "Q_bag_G_mc column (currently all NA)",SMLL),
     pl("Bug identified, not yet re-run",SMLL)],
    [pri("Medium",AMBER),
     pl("GDP legend fix: one-line pt.pch removal to regenerate combined figure",SMLL),
     pl("Trivial — pending",SMLL)],
    [pri("Medium",AMBER),
     pl("Integrate addendum (Section 13) into main proposal: "
        "decide placement (Sec 7 extension vs standalone appendix)",SMLL),
     pl("Decision pending",SMLL)],
    [pri("Medium",AMBER),
     pl("Practitioner rule of thumb: formalise k* &#8712; "
        "[T<super>1/3</super>, T<super>1/2</super>] as a documented recommendation",SMLL),
     pl("Needs paper section",SMLL)],
    [pri("Low",DGREY),
     pl("Fix GDP script bugs for final results: "
        "exclude=1845, ak[&#8722;c(1:2),], dm.test NAs",SMLL),
     pl("Identified, deferred",SMLL)],
]
open_tbl = Table(open_rows, colWidths=[0.65*inch, 4.05*inch, 2.05*inch], repeatRows=1)
open_tbl.setStyle(TableStyle([
    ("BACKGROUND",   (0,0),(-1,0), BLUE),
    ("FONTNAME",     (0,0),(-1,-1),"Helvetica"),
    ("ALIGN",        (0,0),(-1,-1),"LEFT"),
    ("ALIGN",        (0,0),(0,-1),"CENTER"),
    ("VALIGN",       (0,0),(-1,-1),"TOP"),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[WHITE,LBLUE]),
    ("GRID",         (0,0),(-1,-1), 0.35, MGREY),
    ("TOPPADDING",   (0,0),(-1,-1), 3),
    ("BOTTOMPADDING",(0,0),(-1,-1), 3),
    ("LEFTPADDING",  (0,0),(-1,-1), 4),
    ("RIGHTPADDING", (0,0),(-1,-1), 4),
]))
story.append(open_tbl)

story += [Spacer(1, 0.2*inch), hline(0.5, MGREY),
    Paragraph(
        "Generated 4 June 2026 &nbsp;|&nbsp; "
        "Exogenous PID 71717 (joint_B_H_grid running) &nbsp;|&nbsp; "
        "H=6 scaling PID active (500 reps, paper profile) &nbsp;|&nbsp; "
        "GDP outputs: GDP_constrained/all_k_results/ &nbsp;|&nbsp; "
        "LaTeX addendum: rsm_powerlaw_section.tex",
        Caption)]

doc.build(story)
print("Done:", OUTFILE)
