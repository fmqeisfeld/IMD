# IMD
/******************************************************************************
*
* IMD -- The ITAP Molecular Dynamics Program
*
* Copyright 1996-2007 Institute for Theoretical and Applied Physics,
* University of Stuttgart, D-70550 Stuttgart
*
******************************************************************************/
// ***********************************************************************************************************
// *               CODE UPDATES
// * ....     : Einige Ablations-Treshold für Al:
// *    0.064-0.07 J/cm^2 (absorbed) aus Experiment, siehe Paper von Zhakovski (entspricht ~ 700 J/m^2)
// *          mit zugehöriger Kratertiefe von 50 nm
// *    200-400 J/m^2 (abs.) aus Arbeit von Dennis nach verschiedenen Autoren
// *    600 J/m^2 (abs.) mit Zhakovski Potential
// *    1200-1900 J/m^2 (abs.) für Spallation mit Zhakovski Potential
// *
// * .....    : TTM-simulationen lassen sich restarten durch Einlesen der ttm-outputs
// *
// * .....    : Explizites lösen der diffusion fuer variables Kappa nach "konservativem" Schema
// *    Diverse matlab-tests haben gezeigt, dass sowohl für den linearn als auch den nichlinearen
// *    fall, die konservative Formulierung einen geringeren Fehler macht als die ausführlichere,
// *    nicht-konservative, d.h. wenn man alle partiellen ableitungen nach der kettenregel
// *          mit aufnimmt. Weiterer Vorteil: Die Konservative Formulierung ist deutlich simpler
// *
// * .....    : filter option, zum löschen von rausfliegenden atomen erweitert. Funktioniert jetzt in alle richtungen
// *    via filter_min_y, filter_max_y, filter_min_x, filter_max_x, etc.
// *
// * .......  : pdecay wurde nun nach imd_integrate.c verschoben. Dadurch spart man sich eine weitere, unnötige schleife
// *    ueber alle atome
// *
// *          : Generelle Code-Umstrukturierungen zwecks Lesbarkeit etc. mit einigen kleinen "helper"-routinen
// *
// * ......   : Maxwell solver fuer 2D und 1D steht jetzt und funkioniert auch ->Siehe imd_fdtd.c
// *
// * 12.08.18: Dichte wird aus zahl der nachbaratomen pro cutoff-radius-spähre berechnet --> deutlich homogener
// *         weil es dadurch keine rolle mehr spielt ob die Nachbarzelle evtl. eine dicht-gepackte ebene mehr
// *         oder weniger hat. Dazu musste in die struct md_cell eine weitere variable numneighs aufgenommen werden
// *         die auch noch via send-cells bei der kommunikation nicht vergessen werden darf (an vielen stellen
// *         in den imd quelldateien)
// *
// * 14.08.18: local-order-parameter eingeführt (LODP): Gibt die "Kristallinität" pro atom an.Daraus wird dann ein
// *         fd-zellen mittelwert berechnet --> erlaubt die verwendung von separaten wide-range functionen
// *         fuer amorphe und kristalline phase (bisher nur fcc). ACHTUNG: hierfur musste das standardmäßig
// *         aktivierte NEWTON deaktiviert werden (also actio=reactio in Kraf-schleife).
// *         Dadurch ist der code natuerlich etwas langsamer. Außerdem ist LODP bisher nur für Zellenlisten OHNE
// *         NBL-Option eingebaut.
// *
// * 01.09.18: electron-ionen energie-kopplung modifiziert (in do_DIFF): xi wird nun nicht auf Atomzahl der
// *         FD-Zelle, sondern auf "lokale Umgebungsdichte" bezogen. Dies ist in der Routine genauer erläutert
// *
// * 06.09.18: Dirichlet Randbedinungen fuer die elek-temp in 2D an den Raendern der Probe
// *         Achtung:koennte zu stark kuehlen --> großer temp.gradient --> noch schnellere kuehlung
// *         d.h. : Drauf achten, dass die Probe groß genug ist
// *         Noch besser wäre natürlich das TTM-Gitter (über MD-Domäne hinaus) zu erweitern, sodass
// *         die Temperatur "natürlich" diffundieren kann
// *         Neuer Parameter: dirichlet_surfx  gibt die x-Position der Oberfläche an.
// *              sodass die RB nur für Zellen mit x>dirichlet_surfx aktiv werden
// *        Dadurch wird verhindert, dass das Ablatierte Material gekühlt wird
// *
// * 10.10.18: automaische berechnung des maximalen stabilen timesteps fuer diffusion (&advection)
// *         nach von Neumann. ACHTUNG: streng genommen gilt Fourier-analyse nur fuer lineare dgl
// *         wird aber trotzdem auch fuer nichtlin. dgl herangezogen
// *
// * 17.11.18: Mehrere Routinen fuer Interpolation (bilinear,bikubisch, beliebig quadrilateral)
// *         Zweck: Wide-range funktionen aus tabellen statt on-the-fly computation --> Kein Hardcoding fuer
// *         neues Material noetig (Auch die funktioniellen Formen können bei andrem Material anders aussehen, d.h.
// *         es gibt keine standart-funktion der man nur einen satz von fitting-params geben muss,..leider)
// *         --> Größere Flexibilität + man vermeidet teure funktionen wie sin,cos,exp, etc..
// *
// *	
// *
// * 27.11.18: imd_mpiio.c: Input und output via mpiio-schnittstelle. Für mpiio-input -->parallel_input=2
// *         Für mpiio-output: --->parallel_output=3 oder 4
// *         Es wird nun außerdem die Simulationsbox mit rausgeschrieben. Beim kompilieren die option MPIIO mit angeben
// *         Nachtrag: Auch Fortsetzen von Simulationen aus binären *.mpiio-checkpoint-files funzt nun.
// *         Nice to have: Für die Möglichkeit erweitern dass man sich mit single-precision zufrieden gibt
// *
// *
// * 02.12.18: Advection-"Solver" implementiert (ADVMODE=2), welcher den Atomaren Teilchenstrom benutzt (Bisher nur 1D und 2D!!!).
// *         Dabei wird in jedem md-step aus den aktuellen Geschwindigkeiten und Beschleunigungen der Atome vorausgesagt,
// *         wieviele Atome ihre TTM-Zelle verlassen und wohin sie "fließen". Bei dieser Methode ist die inter-prozessor-
// *         kommunikation deutlich aufwändiger als bei der vorherigen, da hier auch Kanten/Ecken-Nachbarn kommunizieren
// *         müssen (Nicht nur an den Seiten der "CPU-Würfel", sondern auch über die Kanten bzw. in 3D auch über die Ecken)
// *         Im Prinzip ist diese Methode dafür genauer als über die Schwerpunkts-Geschw.Vektoren zu gehen
// *         TODO: Kommunikation effizienter machen! Momentan wird komplettes flux-array kommuniziert (in 2D sind das 8 Einträge)
// *         Obwohl nicht alle noetig sind
// *
// * 04.12.18: Paramter shiftx_front und shiftx_rear:in Angstrom
// *         Dadurch kann vor bzw. hinter der Probe Vakuum eingebaut werden. Bei großen Proben ist dies praktischer als
// *         das per z.B. AWK zu machen, da es unter Umständen sehr lange dauern kann. IMD erledigt das ganze
// *         parallelisiert. U
// *         UPDATE: Nun funktioniert das auch direkt in imd_generate.c, also während die Probe erstellt wird.
// *         UPDATE: shifty_front und shifty_rear hinzugefügt
// *
// * 06.12.18: Umgebungsdichte-Berechnung auch für imd_forces_nbl.c eingebaut damit neighlist genutzt werden kann -->SPEEDUP!!
// *         Achtung: Damit die Filter-Option auch mit Nachbarlisten funzt, muss nach jedem loeschen have_valid_nbl=0 gesetzt werden,
// *         sodass die neigh-liste im nächsten md-step aktualisiert wird (geschieht automatisch).
// *         Außerdem ist eine vektorisierung der Kraftschleife einfacher wenn man über neigh-list geht
// *         Nice to have: Die Doppelschleife in eine einzelne Schleife umschreiben --> noch effizientere Vektorisierung möglich
// *         Dazu einfach aus der Neigh-List eine Pair-list bauen (jedes Paar darf nur 1x auftauchen)
// *
// * 14.12.18: Trikubische Interpolation eingebaut. Praktisch für elektronische Parameter, die
// *         abhängen von Dichte, Te und Ti  z.B. kappa,gamma, epsilon, effektive Elek-Koll.freq. etc.
// *         Siehe dazu die Routinen in imd_interpol.c mit dem Prefix "tricub_"
// *         Sparsam verwenden, da sehr Speicher & Rechenintensiv (langsam) trotz guter Vektorisierbarkeit.
// *         (Hier findet eine Multiplikation mit einer 64x64 Matrix statt)
// *         Wenn möglich lieber Hardcoding.
// *
// * 19.12.18: Maxwell-Solver verbessert.
// *         Die vom E-Feld dissipierte Leistungsdichte wird nun nicht mehr harmonisch genähert
// *         (d.h. mittels Slowly-Varying Envelope Approx.- Diese Näherung steckt übrigens auch im Lambert-Beer'schen Gesetz)
// *         Sondern direkt berechnet
// *         D.h. man erhält die "instantane" Leistungsdichte statt der zuvor, über mehrere optsche Zyklen gemittelten Leistungsdichte.
// *         UPDATE: zusätzlich zum reinen Drude-Medium gibt es nun auch die Option für ein Drude-Lorentz Medium
// *		 UPDATE: Reines Drude-Medium entfernt, ergibt sich jedoch als Grenzfall des Lorentz-Modells mit 2 Oszillatoren, bei dem
// *				 eine der Resonanzfrequenzen = 0. 
// *
// * 20.12.18: MPI-output nun auch für ttm-files, da bei großen Simulationen sonst sehr viel kommuniziert werden muss,
// *         was entsprechend lange dauert.
// *         Wichtig hierbei ist, dass die procs ihren Zellen einfach hintereinander
// *         rausschreiben, also nicht nach Ort sortiert. Deswegen hat gnuplot Probleme hieraus eine heatmap zu plotten.
// *         Die Zeilen müssen zuvor umsortiert werden mittels
// *
// *         sort -g -k1,1 -k2,2 -k3,3 laser.1.ttm > laser.1.sorted.ttm
// *
// * 27.12.18: Instantane Absorbierte Leistungsdichte korrigiert (Siehe imd_fdtd.c)
// *         ACHTUNG: bisher ist 1D-bzw. 2D-Option HARDCODING mittels #define FDTD1D bzw. #define FDTD2D in imd_fdtd.c
// *
// * 31.01.19: Randbedingungen fuer Godunov-Flux-solver in y/z nichtperiodisch korrigiert. Hierbei wird angenommen,
// *         dass außerhalb der Simulationsdomäne die Temperaturen = 0 sind (nicht optimal).
// *         Besser wären absorbierende RB ähnlich wie im Maxwell-Solver
// *		 ACHTUNG: Dieser Solver wurde entfernt..Advektiert wird nur noch über die erste Variante, bei dem aus/einsträmende 
// *		 Teilchen direkt gezählt werden
// *		
// *
// *04.06.19:  Transfer-Matrix methode hinzugefügt (helmholtz-solver). Siehe imd_tmm.c. Nur 1D-Simulationen möglich!
// *         Berechnung der Felderverteilung ist NICHT parallisiert. Sollte bei Quasi-1D Simulationen aber keine Rolle spielen.
// *
// *22.07.19:  Schockabsorbierende Randbedingungen nach Comput. Mech. 50:645-655 (2012) für 100-orienterten fcc-kristall
// *         in x und y-Richtung implementiert. Siehe imd_nrb.c
// *         TODO: Randbedinungen "fortsetzbar" machen, d.h. zu Beginn (nachdem die Nachbaratome der Randatome identifiziert sind)
// *               alle nötigen Infos wie indizes der Nachbaratome sowie "Gleichgewichtspositionen" rausschreiben, zum späteren Einlesen.
// *         UPDATE: erledigt. Details hierzu in imd_nrb.c
// *
// *02.08.19:  FILTER-Option so erweitert, dass nun vor dem Löschen der Atome außerhelb der Filter-Grenze geprüft wird,
// *         ob das Atom "irgendwie" mit einem Atom verbunden ist, dessen Position noch innerhalb der Grenzen liegt (Cluster-Analyse).
// *         D.h. es wird für jedes vermeintlich zu löschende Atom die Nachbarkette geprüft, bis ein Atom gefunden wird,
// *         das nicht gelöscht werden darf. Algorithmisch ist die mittels eines rekursiven Schemas umgesetzt, was auch bedeutet
// *         dass es je nach verfügbarem Speicher eine bestimmte obere Grenze für die Rekursionstiefe gesetzt werden muss
// *         (Nicht zuletzt wegen der Performance)
// *         Diese Grenze wurde auf 1000 gesetzt. In der Regel ist diese Rekursionstiefe nicht nötig, denn sobald
// *         ein Nachbarpfad in eine "Sackgasse" läuft, d.h. dass es in der Umgebung keine Nachbaratome gibt, die nicht
// *         bereits vorher besucht worden sind, wird die Suche beim vorherigen Atom fortgesetzt usw.
// *         Sobald die Kette auf diese Weise wieder am Ursprungsatom angelangt ist, kann man sicher sein, dass es
// *         tatsächlich keine Verbindung zu einem Atom gibt, dass nicht gelöscht werden darf.
// *         Einzelheiten in imd_filter.c
// *
// * TODO:
// *       Etwas stimmt nicht mit umgebungsdichte-berechnung....alle paar steps abrupte änderung...
// *        --> I.was in Kommunikation etc.
// ***********************************************************************************************************