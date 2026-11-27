# CIR Numerical Precision Lab – gemeinsamer Projektprompt für Claude Code

Wir starten ein neues Projekt komplett von Grund auf. Das aktuelle Verzeichnis ist leer.

Ich möchte mit dir eine kleine Anwendung entwickeln, die untersucht, welchen Einfluss die numerische Rechengenauigkeit auf die Monte-Carlo-Bewertung einer bermudischen Call-Option auf einen Coupon-Bond innerhalb eines Cox-Ingersoll-Ross-Short-Rate-Modells hat.

Konkret sollen zwei ansonsten möglichst identische Implementierungen miteinander verglichen werden:

1. Berechnung mit `double`
2. Berechnung mit `BigDecimal`

Das primäre Ziel besteht nicht nur darin, die resultierenden Optionspreise zu vergleichen, sondern insbesondere darin zu untersuchen, ob und an welchen Stellen die höhere arithmetische Präzision von `BigDecimal` zu messbaren Unterschieden führt.

Ein besonderer Schwerpunkt soll auf der Rückwärtsinduktion und der Schätzung der bedingten Fortsetzungswerte liegen. Als primäre Methode sollen hierfür zustandsbasierte Quantilcluster verwendet werden. Zusätzlich soll als zweite Methode Longstaff-Schwartz Least-Squares Monte Carlo vorgesehen werden, damit später untersucht werden kann, ob die Wahl des Verfahrens zur Schätzung der Fortsetzungswerte einen größeren Einfluss besitzt als die numerische Repräsentation durch `double` beziehungsweise `BigDecimal`.

Handle als erfahrener Quantitative Developer, numerischer Mathematiker und Java Software Engineer. Das Projekt ist ein wissenschaftlicher Prototyp, keine produktive Bewertungsbibliothek. Arbeite ergebnisoffen, nachvollziehbar und reproduzierbar.

**Priorität:** Das Hauptziel ist die bermudische Call-Option auf eine Kuponanleihe mit Quantilclustern, einschließlich der numerischen Stabilität der Rückwärtsinduktion. Europäische Optionen und analytische Bewertungsverfahren dienen als Kontroll- und Validierungsexperimente. Longstaff-Schwartz, Portfolios und vertiefte Performance-Experimente folgen schrittweise. Vermeide, dass optionale Erweiterungen den ersten fachlich nutzbaren Prototyp verzögern.

Untersuche neben der Genauigkeit insbesondere die Berechnungsdauer und den Ressourcenbedarf, zunächst bei einzelnen Verträgen und später bei ganzen Vertragsbeständen. Die Anwendung soll eine Entscheidung darüber ermöglichen, ob und an welchen Stellen sich höhere Präzision gegenüber `double` praktisch lohnt.

## 1. CIR-Modell

Unter dem risikoneutralen Maß $Q$ soll die Short Rate $r$ der SDE

$$
dr_t=\theta(\mu-r_t)\,dt+\sigma\sqrt{r_t}\,dW_t
$$

folgen, mit

$$
\theta>0,\qquad \mu>0,\qquad \sigma>0,\qquad r_0>0.
$$

Alle vier Größen sollen als Modellparameter übergeben werden können.

Im Standardmodus müssen die Parameter die Feller-Bedingung

$$
2\theta\mu\ge\sigma^2
$$

erfüllen. Andernfalls sollen die Modellparameter im Standardmodus bereits bei ihrer Erzeugung mit einer geeigneten Exception abgelehnt werden.

Für spätere Robustheitsexperimente soll ein expliziter Modus `ALLOW_VIOLATION` vorgesehen werden, während `STRICT` die Standardeinstellung bleibt. Die übrigen Parameterrestriktionen gelten in beiden Modi. Eine Verletzung der Feller-Bedingung bedeutet nicht, dass der kontinuierliche CIR-Prozess negative Werte erzeugt; sie ist von negativen Werten der Euler-Approximation zu unterscheiden.

Verwende die hier eingeführte Notation (theta = mean reversion speed, mu = long-term mean) durchgängig. In Java sollen sprechende Namen wie `meanReversionSpeed`, `longTermMean`, `volatility` und `initialRate` verwendet werden.

## 2. Underlying und bermudische Option

Das Underlying ist ein Strom zukünftiger positiver Cashflows

$$
C_1,\ldots,C_n
$$

zu den Zeitpunkten

$$
0<t_1<\ldots<t_n.
$$

Insbesondere soll unkompliziert eine einfache Kuponanleihe mit regulären jährlichen Kupons erzeugt werden können:

$$
C_1=\ldots=C_{n-1}=cN,\qquad C_n=(1+c)N,
$$

mit Nominal $N$ und Kuponzinssatz $c$.

Die Option ist eine bermudische Call-Option mit Strike $K>0$ und einer vorgegebenen Menge von Exercise Dates

$$
0<e_1<e_2<\ldots<e_m<t_n.
$$

An jedem Exercise Date kann der dann noch verbleibende Cashflow-Strom zum Strike $K$ erworben werden.

Der Wert des Underlyings zum Exercise Date $e_j$ ist

$$
B_{e_j}
=
\sum_{t_i>e_j} C_iP(e_j,t_i).
$$

Der Immediate Exercise Value ist damit

$$
E_j=(B_{e_j}-K)^+.
$$

Die Zerobondpreise $P(s,T)$ sollen aus der analytischen CIR-Zerobondpreisformel ermittelt werden. Die Zerobondpreise selbst sollen nicht durch zusätzliche Monte-Carlo-Simulation approximiert werden.

Der Immediate Exercise Value ist ausdrücklich vom Continuation Value der Option zu unterscheiden. Die Notation $C_i$ für Bond-Cashflows und $C_j$ für den Continuation Value ist aus dem jeweiligen Kontext zu unterscheiden; im Code eindeutige Namen wie `bondCashflow` und `continuationValue` verwenden.

Am Exercise Date $e_j$ ist der theoretische Continuation Value der bedingte Erwartungswert des auf $e_j$ diskontierten zukünftigen Optionscashflows unter der ab diesem Zeitpunkt optimalen Exercise Policy:

$$
C_j
=
E^Q\left[
Y_j
\mid
\mathcal F_{e_j}
\right].
$$

Dabei bezeichnet $Y_j$ den auf $e_j$ diskontierten Cashflow der nach aktuellem Stand der Rückwärtsinduktion optimalen späteren Ausübung.

Die relevante Bellman-Entscheidung lautet somit

$$
V_j=\max(E_j,C_j).
$$

Da $C_j$ nicht analytisch bestimmt werden soll, wird er aus den simulierten Pfaden geschätzt.

Für den ersten Prototyp genügen Zeitangaben in Jahren und ein Zeitraster, das sämtliche Exercise Dates exakt enthält. Kalenderlogik und Day-Count-Conventions sind nicht erforderlich. Die obige Cashflow-Konvention ist verbindlich: Zahlungen am Exercise Date selbst gehören nicht zum verbleibenden Bond. Der Strike ist ein Geldbetrag in derselben Einheit wie der Bondwert (nicht automatisch ein prozentualer Strike). Vermeide eine doppelte Nominalskalierung.

Ein allgemeines Cashflow-Modell soll auch einen einzelnen Cashflow und ein einzelnes Exercise Date unterstützen. So lassen sich die europäischen Referenzfälle ohne separates Vertragsmodell darstellen.

## 3. Rückwärtsinduktion

Die optimale Exercise Policy soll durch Rückwärtsinduktion über die Exercise Dates bestimmt werden.

Am letzten Exercise Date $e_m$ gibt es keine spätere Exercise Opportunity mehr. Daher gilt

$$
V_m=E_m.
$$

An jedem vorherigen Exercise Date $e_j$ soll für jeden relevanten Pfad der bereits aus der Rückwärtsinduktion bekannte zukünftige optimale Cashflow auf $e_j$ diskontiert werden.

Für Pfad $p$ sei $\tau^{(p)}>e_j$ der nach aktuellem Stand der Rückwärtsinduktion gewählte spätere Exercise-Zeitpunkt. Dann ist das Target für die Schätzung des Fortsetzungswerts

$$
Y_j^{(p)}
=
D^{(p)}(e_j,\tau^{(p)})
E_{\tau^{(p)}}^{(p)}.
$$

Aus diesen Targets wird der bedingte Continuation Value

$$
\widehat C_j
$$

geschätzt.

Anschließend erfolgt pfadweise die Entscheidung

$$
E_j^{(p)}>\widehat C_j^{(p)}
\quad\Longrightarrow\quad
\text{Exercise}.
$$

Wird auf Pfad $p$ zu $e_j$ ausgeübt, ersetzt der Exercise Cashflow zum Zeitpunkt $e_j$ den bislang vorgesehenen späteren Cashflow dieses Pfades.

Dadurch beeinflusst eine Exercise Decision an einem späteren Exercise Date gegebenenfalls die Schätzung der Fortsetzungswerte und Exercise Decisions an früheren Exercise Dates.

Diese Rückwärtsinduktion soll in der Implementierung möglichst transparent nachvollziehbar sein.

Definiere explizit den Fall, dass auf einem Pfad keine spätere Ausübung erfolgt: Der zukünftige Optionscashflow und sein diskontiertes Target sind dann null. Verwende bei Gleichstand eine deterministische Konvention (Vorschlag: Continue); übe nur bei positivem Immediate Exercise Value aus. Unterscheide die während der Rückwärtsinduktion hypothetisch bestimmten Entscheidungen von der in Vorwärtsrichtung tatsächlich realisierten einmaligen Ausübung.

Speichere für die Rückwärtsinduktion vorzugsweise den jeweils ausgewählten späteren Ausübungszeitpunkt und Cashflow je Pfad, nicht einen unkontrolliert wachsenden Baum von Zwischenwerten.

## 4. Primäre Schätzung der Fortsetzungswerte: Quantilcluster

Als primäre Methode zur Schätzung der bedingten Fortsetzungswerte sollen zustandsbasierte Quantilcluster verwendet werden.

Da das CIR-Modell eindimensional und Markovsch ist, soll zunächst die aktuelle Short Rate

$$
r_{e_j}
$$

als Zustandsvariable für die Clusterbildung verwendet werden.

Im Entwicklungsplan soll zusätzlich kurz beurteilt werden, ob alternativ der aktuelle Bondwert $B_{e_j}$ eine sinnvollere Clustervariable wäre. Für den ersten Prototypen soll jedoch nur eine Zustandsvariable verwendet werden.

An jedem Exercise Date sollen die relevanten Pfade anhand ihrer Zustandsvariable sortiert und anschließend in möglichst gleich große Quantilcluster aufgeteilt werden.

Beispielsweise können bei hinreichend vielen Pfaden Größenordnungen von 20 bis 50 Clustern untersucht werden. Die konkrete Wahl soll konfigurierbar sein und im Entwicklungsplan begründet werden.

Für einen Cluster $G_k$ soll der Fortsetzungswert als Mittelwert der diskontierten zukünftigen Cashflows der diesem Cluster zugeordneten Pfade geschätzt werden:

$$
\widehat C_{j,k}
=
\frac{1}{|G_k|}
\sum_{p\in G_k}Y_j^{(p)}.
$$

Für einen Pfad $p\in G_k$ gilt damit

$$
\widehat C_j^{(p)}=\widehat C_{j,k}.
$$

Der Fortsetzungswert ist bei dieser Methode also eine lokale, stückweise konstante Schätzung des bedingten Erwartungswerts.

Die Anzahl der Cluster soll konfigurierbar sein.

Es soll außerdem berücksichtigt werden, dass eine zu geringe Clusterzahl zu einem erhöhten Approximationsbias führt, während eine zu hohe Clusterzahl die statistische Varianz der lokalen Mittelwerte erhöht. Leere oder zu kleine Cluster sollen durch die Verwendung von Quantilclustern möglichst vermieden werden.

Im Entwicklungsplan soll eine einfache Mindestanforderung an die Anzahl der Pfade pro Cluster vorgeschlagen werden.

## 5. Relevante Pfade für die Schätzung des Fortsetzungswerts

Im Entwicklungsplan soll ausdrücklich festgelegt und begründet werden, welche Pfade zur Schätzung des Fortsetzungswerts herangezogen werden.

Insbesondere soll diskutiert werden, ob analog zu einer typischen Longstaff-Schwartz-Implementierung nur aktuell In-the-Money-Pfade mit

$$
E_j^{(p)}>0
$$

verwendet werden oder ob für die Quantilcluster alle noch lebenden Pfade sinnvoller sind.

Diese Entscheidung soll nicht implizit getroffen werden, da sie einen relevanten Einfluss auf die lokale Schätzung des Continuation Value haben kann.

Für den ersten Prototypen soll eine einfache und einheitliche Konvention gewählt und durch Tests dokumentiert werden.

Als Ausgangspunkt ist eine Schätzung auf allen Trainingspfaden am betreffenden Exercise Date zu prüfen; die Ausübung selbst bleibt auf positive Immediate Exercise Values beschränkt. Bei der Rückwärtsinduktion darf die Pfadauswahl nicht von einer erst später berechneten früheren Ausübungsentscheidung abhängen. Dokumentiere alternativ die ITM-Auswahl und ihre Auswirkungen. Für den kontrollierten Vergleich müssen auch die zur Schätzung zugelassenen Pfadindizes identisch sein, oder Unterschiede in der Auswahl als eigenständiger Effekt ausgewiesen werden.

## 6. Gemeinsame und getrennte Quantilcluster

Der `double`-gegen-`BigDecimal`-Vergleich soll zwei unterschiedliche Modi unterstützen.

### Modus A: gemeinsame Quantilcluster

Im kontrollierten Vergleich sollen beide Implementierungen dieselbe Clusterzuordnung der Pfade verwenden.

Die Cluster sollen also dieselben Pfadindizes enthalten. Innerhalb dieser identischen Cluster werden anschließend die jeweiligen diskontierten zukünftigen Cashflows mit `double` beziehungsweise `BigDecimal` verarbeitet.

Dieser Modus soll möglichst stark isolieren, welche Auswirkungen die unterschiedliche Arithmetik auf

- Pfadsimulation,
- Diskontierung,
- Summation,
- Mittelwertbildung,
- Continuation Values und
- Exercise Decisions

hat.

### Modus B: getrennte Quantilcluster

In einem zweiten End-to-End-Vergleich soll jede Implementierung ihre Quantilcluster selbst anhand ihrer jeweils berechneten Zustandswerte bilden.

Damit dürfen geringfügige Unterschiede in den simulierten Short Rates auch zu Unterschieden in

- Sortierreihenfolge,
- Clusterzuordnung,
- Cluster-Mittelwerten,
- Continuation Values und
- Exercise Decisions

führen.

Dieser Modus soll zeigen, ob kleine numerische Abweichungen durch die diskreten Entscheidungen der Clusterbildung und der Rückwärtsinduktion verstärkt werden.

Die Anwendung soll beide Modi sauber voneinander unterscheiden und ihre Ergebnisse getrennt ausweisen.

Lege für Modus A fest, mit welcher Zustandsdarstellung die gemeinsamen Cluster erzeugt werden (z. B. die Double-Zustände als feste Referenz). Verwende eine deterministische Sortierung mit Pfad-ID als Tie-Breaker und dokumentiere die Behandlung identischer Zustandswerte, insbesondere an der Grenze zweier Cluster. Definiere Mindestclustergröße, Behandlung kleiner Stichproben sowie die Intervallgrenzen eindeutig.

Gemeinsame Trainingscluster bedeuten noch nicht identische Exercise Policies. Optional soll später ein zusätzlicher Kontrolllauf mit eingefrorener gemeinsamer Exercise Policy den Einfluss der reinen Cashflow-Arithmetik von dem Einfluss abweichender Entscheidungen trennen. Für unabhängige Bewertungspfade ist im kontrollierten Modus auch die Zuordnung zu den trainierten Clustern konsistent zu halten.

## 7. Optionale zweite Methode: Longstaff-Schwartz

Zusätzlich zur primären Cluster-Methode soll die Architektur so geplant werden, dass als zweite Methode Longstaff-Schwartz Least-Squares Monte Carlo verwendet werden kann.

Für den ersten vollständigen Prototypen darf Quantil-Clustering jedoch priorisiert und Longstaff-Schwartz anschließend in einem getrennten Entwicklungsschritt ergänzt werden.

Beim Longstaff-Schwartz-Verfahren soll derselbe zukünftige diskontierte optimale Cashflow

$$
Y_j^{(p)}
$$

als Regression Target verwendet werden.

Als erste einfache Regressionsbasis kann beispielsweise

$$
1,\quad r,\quad r^2
$$

verwendet werden.

Alternativ soll im Plan kurz diskutiert werden, ob

$$
1,\quad B,\quad B^2
$$

mit dem aktuellen Bondwert als erklärender Variable sinnvoller wäre.

Es sollen zunächst keine komplexen Machine-Learning-Verfahren oder unnötig großen Regressionsbasen verwendet werden.

Für die Lösung des Regressionsproblems soll keine explizite Matrixinversion verwendet werden. Im Entwicklungsplan soll insbesondere eine QR-basierte Lösung als mögliche robuste Standardmethode betrachtet werden.

Ein späterer Vergleich soll insbesondere die Größen

$$
\Delta_{\mathrm{precision}}
=
V_{\mathrm{BigDecimal}}-V_{\mathrm{double}}
$$

innerhalb derselben Fortsetzungswertmethode und

$$
\Delta_{\mathrm{method}}
=
V_{\mathrm{LSM}}-V_{\mathrm{Cluster}}
$$

bei derselben numerischen Repräsentation untersuchen.

Damit soll geprüft werden, ob die Wahl der Methode zur Schätzung des Fortsetzungswerts einen größeren Einfluss hat als die Wahl zwischen `double` und `BigDecimal`.

## 8. Monte-Carlo-Simulation

Die Short Rate soll mit einem Euler-Maruyama-Schema simuliert werden.

Die beiden Implementierungen sollen hinsichtlich

- Modell,
- Diskretisierung,
- Anzahl der Pfade,
- Zeitraster,
- Exercise Dates,
- Cashflows,
- Seed,
- standardnormalverteilten Zufallszahlen und
- Rückwärtsinduktionslogik

so identisch wie möglich sein.

Der einzige beabsichtigte Unterschied soll die verwendete numerische Darstellung und Arithmetik sein.

### Variante A: `double`

Sämtliche numerischen Größen der eigentlichen Simulation und Bewertung sollen als `double` dargestellt werden.

### Variante B: `BigDecimal`

Modellparameter, Zustandsgrößen, Cashflows, Diskontfaktoren sowie die arithmetischen Operationen der Simulation und Bewertung sollen soweit sinnvoll mit `BigDecimal` umgesetzt werden.

Im Entwicklungsplan soll ausdrücklich untersucht werden, an welchen Stellen eine vollständig auf `BigDecimal` basierende Implementierung problematisch ist, insbesondere bei

- Quadratwurzeln,
- Exponentialfunktionen,
- Erzeugung standardnormalverteilter Zufallszahlen,
- Sortierung und Clustergrenzen sowie
- gegebenenfalls später der linearen Algebra für Longstaff-Schwartz.

Falls Konvertierungen zwischen `double` und `BigDecimal` sinnvoll oder notwendig sind, sollen diese auf möglichst kleine und klar abgegrenzte Stellen beschränkt und dokumentiert werden.

Für `BigDecimal` soll eine explizite und konfigurierbare `MathContext`-Strategie vorgesehen werden. Für den ersten Prototypen darf beispielsweise `MathContext.DECIMAL128` als Standard betrachtet werden.

### Konkrete Präzisionsstrategie

Verwende als primäre Variante `MathContext.DECIMAL128` (34 signifikante Stellen, `HALF_EVEN`). Präzision und Rundungsmodus sollen konfigurierbar sein. Verwende für ausgewählte Referenzläufe 64 Stellen und überprüfe deren Stabilität mit 100 Stellen oder höherer Präzision, soweit dies praktikabel ist.

Für mathematische Funktionen ist die Bibliothek `ch.obermuhlner:big-math` als Ausgangspunkt vorgesehen. Prüfe eine verfügbare und kompatible Version. Nutze insbesondere hochpräzise Exponentialfunktionen, Quadratwurzeln und Logarithmen. Implementiere elementare mathematische Funktionen nicht unnötig selbst. Vermeide in der eigentlichen BigDecimal-Bewertung eine Rückkonvertierung nach `double` (z. B. `Math.exp(value.doubleValue())`). Ausnahmen müssen explizit abgegrenzt, begründet und als solche ausgewiesen werden.

- Fachliche Dezimaleingaben aus Strings einlesen und beide Varianten aus derselben Eingabekonfiguration ableiten.
- `MathContext` bei den relevanten Rechenoperationen konsistent anwenden.
- Arbeitspräzision, interne Schutzstellen bei Funktionsapproximationen und Ausgaberundung unterscheiden.
- Keine vorzeitige Rundung auf feste Nachkommastellen während der Bewertung.
- `BigDecimal`-Differenzen ohne vorzeitige Umwandlung beider Ergebnisse nach `double` berechnen.
- Eine höherpräzise Rechnung ist eine numerische Referenz desselben Algorithmus, kein Beweis für die Richtigkeit der Bewertungsmethode.

Spätere hybride Kontrollvarianten können `double` mit Kahan- oder paarweiser Summation sowie `double`-Simulation mit `BigDecimal`-Aggregation umfassen. Sie dienen der Lokalisierung von Fehlerquellen, nicht als Ersatz für die primäre durchgängige BigDecimal-Variante.

## 9. Negative Werte beim Euler-Maruyama-Schema

Ein fachlich wichtiger Punkt ist, dass das naive Euler-Maruyama-Verfahren trotz erfüllter Feller-Bedingung diskret negative Short Rates erzeugen kann.

Dadurch wird insbesondere

$$
\sqrt{r_t}
$$

problematisch.

Dieses Problem soll im Entwicklungsplan ausdrücklich analysiert werden.

Für den ersten Prototypen soll eine einfache und transparente Behandlung vorgeschlagen werden. Es soll nicht stillschweigend auf ein anderes Diskretisierungsschema gewechselt werden.

Insbesondere soll deutlich zwischen Eigenschaften des kontinuierlichen CIR-Prozesses und Artefakten seiner zeitdiskreten Approximation unterschieden werden.

Als bevorzugten Ausgangspunkt prüfe Full-Truncation Euler. Definiere mit $x^+=\max(x,0)$ explizit:

$$
x_{k+1}=x_k+\theta(\mu-x_k^+)\Delta t+\sigma\sqrt{x_k^+}\,\Delta W_k.
$$

Unterscheide den intern fortgeschriebenen Zustand $x_k$ vom nichtnegativen effektiven Zinssatz $x_k^+$. Dokumentiere, welcher Zustand für Bondbewertung, Clustering und Diskontierung verwendet wird (Vorschlag: der effektive Zinssatz). Eine Projektion des gespeicherten Zustands nach jedem Schritt wäre eine andere Variante und muss entsprechend gekennzeichnet werden. Zähle negative Zwischenwerte und die Anwendung der Korrektur in beiden numerischen Varianten.

Milstein und exakte CIR-Übergangssimulation sind spätere, nicht verpflichtende Erweiterungen. Exakte Übergangswerte an diskreten Zeitpunkten sind nicht gleichbedeutend mit einer exakten Simulation des integrierten Zinspfades.

## 10. Diskontierung entlang der Monte-Carlo-Pfade

Die pfadabhängigen Diskontfaktoren zwischen zwei Zeitpunkten sollen aus der simulierten Short Rate bestimmt werden.

Das Integral

$$
\int_a^b r_s\,ds
$$

muss aufgrund der zeitdiskreten Simulation numerisch approximiert werden.

Beide Implementierungen müssen exakt dieselbe Integrationsregel und dasselbe Zeitraster verwenden.

Im Entwicklungsplan soll eine einfache Integrationsregel für den Prototypen vorgeschlagen und begründet werden.

Beginne mit der Linksrechteckregel auf dem gemeinsamen Zeitraster, sofern keine fachlichen Gründe dagegen sprechen. Die Trapezregel ist eine spätere Kontrollvariante. Speichere bei Bedarf kumulierte Zinsintegrale an Exercise Dates, aus deren Differenzen die Diskontfaktoren zwischen zwei Exercise Dates gebildet werden können. Beide Rechenvarianten müssen dieselbe Regel und denselben effektiven Zustandsbegriff verwenden.

## 11. Reproduzierbarkeit und Common Random Numbers

Der Vergleich zwischen `double` und `BigDecimal` soll möglichst wenig durch gewöhnliches Monte-Carlo-Rauschen beeinflusst werden.

Deshalb sollen beide Varianten exakt dieselben standardnormalverteilten Zufallszahlen verwenden.

Die Zufallszahlenerzeugung soll konzeptionell von der Pfadsimulation getrennt werden.

Für einen festen Seed müssen die verwendeten Zufallszahlen und damit die Simulation reproduzierbar sein.

Ziel ist ein Experiment mit Common Random Numbers, bei dem Unterschiede zwischen beiden Implementierungen möglichst direkt auf ihre numerische Arithmetik zurückgeführt werden können.

Verwende die gleichen tatsächlich erzeugten Normalvariablen (nicht bloß gleiche Seeds). Eine gemeinsame Double-basierte Zufallsquelle ist für das Hauptexperiment zulässig; die Zufallseingaben sind dann auf die Präzision dieser Quelle begrenzt. Für den möglichst isolierten Arithmetikvergleich übertrage die repräsentierten binären Double-Zufallswerte exakt mittels `new BigDecimal(z)` und dokumentiere, wann der MathContext erstmals angewendet wird. `BigDecimal.valueOf(z)` verwendet dagegen die kanonische dezimale Stringdarstellung und ist eine gesonderte Konvertierungskonvention. Dies ist von der String-basierten Einlesung fachlicher Modellparameter zu unterscheiden.

Verwende zunächst einen festen, konkret dokumentierten Generator und Seed. Ergänze später mehrere Seeds zur statistischen Validierung. Bei unterschiedlichen Zeitrastern sollen, soweit möglich, gekoppelte Brownsche Inkremente verwendet werden (grobe Inkremente als Summe feiner Inkremente). Plane speicherschonende Zufallsblöcke und deterministische Pfadindizes ein; bei paralleler Ausführung dürfen Thread-Scheduling und Reihenfolge der Aufrufe die Pfadzuordnung nicht ändern.

## 12. Zielgrößen des Vergleichs

Es sollen nicht nur die finalen Optionspreise verglichen werden.

Mindestens folgende Größen sollen untersucht werden.

### 12.1 Finaler Optionspreis

Für beide Varianten:

- geschätzter Optionspreis,
- Monte-Carlo-Standardfehler,
- Laufzeit.

Zusätzlich:

$$
\Delta_V
=
V_{\mathrm{BigDecimal}}-V_{\mathrm{double}},
$$

sowie absoluter und relativer Preisunterschied.

Der Preisunterschied soll außerdem ins Verhältnis zum Monte-Carlo-Standardfehler gesetzt werden, um beurteilen zu können, ob er gegenüber der statistischen Unsicherheit materiell ist.

### 12.2 Unterschiede der simulierten Short-Rate-Pfade

Da beide Implementierungen dieselben Zufallszahlen verwenden, sollen aggregierte Unterschiede der simulierten Short Rates untersucht werden:

- maximale absolute Abweichung,
- mittlere absolute Abweichung,
- geeignete relative Abweichung,
- Entwicklung dieser Abweichungen mit zunehmender Simulationszeit.

Es sollen hierfür keine unnötig großen Datenmengen gespeichert oder ausgegeben werden.

### 12.3 Unterschiede der Diskontfaktoren

Analog sollen die berechneten Diskontfaktoren verglichen werden:

- maximale absolute Abweichung,
- mittlere absolute Abweichung,
- relative Abweichung.

### 12.4 Unterschiede der Cluster-Fortsetzungswerte

Für jedes Exercise Date und jeden Cluster sollen die mit `double` beziehungsweise `BigDecimal` erhaltenen Schätzungen

$$
\widehat C_{j,k}^{\mathrm{double}}
\quad\text{und}\quad
\widehat C_{j,k}^{\mathrm{BigDecimal}}
$$

verglichen werden.

Mindestens betrachtet werden sollen:

- absolute Differenz,
- relative Differenz,
- Anzahl der Pfade im Cluster,
- Spannweite der Zustandsvariable innerhalb des Clusters,
- Streuung der Regression Targets beziehungsweise zukünftigen diskontierten Cashflows innerhalb des Clusters.

Damit soll erkennbar werden, ob Unterschiede der Fortsetzungswerte tatsächlich aus numerischer Arithmetik stammen oder gegenüber der statistischen Streuung innerhalb eines Clusters vernachlässigbar sind.

### 12.5 Unterschiede der Clusterzuordnung

Im Modus mit getrennten Clustern soll für jedes Exercise Date untersucht werden:

- wie viele Pfade in beiden Varianten demselben Cluster zugeordnet werden,
- wie viele Pfade den Cluster wechseln,
- ob Clusterwechsel überwiegend bei Pfaden nahe einer Quantilgrenze auftreten,
- wie groß die Short-Rate-Abweichung bei solchen Pfaden ist.

Diese Analyse soll zeigen, ob sehr kleine numerische Unterschiede über die diskrete Clusterbildung verstärkt werden.

### 12.6 Unterschiede der Exercise Decisions

Für jedes Exercise Date soll bestimmt werden:

- Anzahl und Anteil der Pfade, auf denen `double` und `BigDecimal` dieselbe Exercise Decision treffen,
- Anzahl und Anteil der Pfade mit unterschiedlicher Entscheidung,
- `double` exercised / `BigDecimal` continued,
- `BigDecimal` exercised / `double` continued.

### 12.7 Abstand zur Exercise Boundary

Für Pfade mit unterschiedlichen Exercise Decisions soll insbesondere

$$
|E_j-\widehat C_j|
$$

untersucht werden.

Dadurch soll geprüft werden, ob abweichende Entscheidungen hauptsächlich auf Pfaden auftreten, die ohnehin sehr nahe an der geschätzten Exercise Boundary liegen.

Als aggregierte Kennzahlen sollen mindestens Median, Mittelwert und Maximum dieses Abstands betrachtet werden.

### 12.8 Rückwirkung durch die Backward Induction

Es soll nach Möglichkeit nachvollziehbar gemacht werden, ob unterschiedliche Exercise Decisions an späteren Exercise Dates zu unterschiedlichen Regression Targets beziehungsweise Cluster-Mittelwerten an früheren Exercise Dates führen.

Ziel ist es zu untersuchen, ob kleine numerische Unterschiede durch die Rückwärtsinduktion abgeschwächt oder verstärkt werden.

### 12.9 Approximative Exercise Boundary

Da das CIR-Modell eindimensional ist, soll untersucht werden, ob für jedes Exercise Date eine approximative Exercise Boundary bezüglich der Short Rate sinnvoll aus den Ergebnissen abgeleitet werden kann.

Die Boundaries der `double`- und `BigDecimal`-Variante sollen miteinander verglichen werden.

Dabei soll berücksichtigt werden, dass die Quantilcluster-Methode einen stückweise konstanten Fortsetzungswert erzeugt, sodass die geschätzte Boundary entsprechend diskret beziehungsweise verrauscht sein kann.

### 12.10 Methodenvergleich Cluster versus Longstaff-Schwartz

Nach Ergänzung der Longstaff-Schwartz-Methode sollen bei identischen Simulationsdaten verglichen werden:

- Optionspreis,
- Standardfehler,
- Laufzeit,
- Exercise Rates je Exercise Date,
- approximative Exercise Boundaries,
- Anzahl beziehungsweise Anteil unterschiedlicher Exercise Decisions.

Besonders interessant ist dabei, ob

$$
|V_{\mathrm{LSM}}-V_{\mathrm{Cluster}}|
$$

größer ist als

$$
|V_{\mathrm{BigDecimal}}-V_{\mathrm{double}}|.
$$

### 12.11 Gepaarte statistische Analyse

Berechne zusätzlich zu den individuellen Standardfehlern und 95%-Konfidenzintervallen die pfadweisen Differenzen der auf den Bewertungszeitpunkt diskontierten Optionscashflows:

$$
d_p=Y_{0,p}^{\mathrm{BigDecimal}}-Y_{0,p}^{\mathrm{double}},\qquad
\overline d=\frac{1}{n}\sum_{p=1}^n d_p,\qquad
SE_d=\frac{s_d}{\sqrt n}.
$$

Berechne ein geeignetes Konfidenzintervall für die gepaarte Preisdifferenz. Bei unabhängig erzeugten Bewertungspfaden und zuvor fixierten Exercise Policies ist dieses bedingt auf die trainierten Policies interpretierbar. Es erfasst nicht automatisch deren Trainingsunsicherheit. Wiederholungen mit unabhängigen Trainings- und Bewertungsseeds dienen der Beurteilung dieser zusätzlichen Unsicherheit.

Berichte absolute und relative Abweichung, Mittelwert, Maximum, RMSE und, wo sinnvoll, ULP-Abstände. Bei Referenzwerten null oder nahe null ist eine relative Abweichung gesondert zu kennzeichnen. Unterscheide numerische, statistische, Zeitdiskretisierungs- und Continuation-Approximationseffekte. Ein Abstand zu einer Referenz enthält möglicherweise mehrere dieser Fehler gleichzeitig.

### 12.12 Qualität der Exercise Policy (Out-of-Sample)

Die erste Implementierung darf die Rückwärtsinduktion auf derselben Pfadmenge demonstrieren. Für belastbare Aussagen sollen die Training Paths zur Bestimmung der Exercise Policy von unabhängigen Valuation Paths zur anschließenden Bewertung getrennt werden. Verwende innerhalb jeder Pfadmenge wiederum Common Random Numbers für beide Arithmetiken. So werden In-Sample-Effekte bzw. Look-ahead-Bias nicht mit numerischen Fehlern verwechselt.

Speichere pro Exercise Date die trainierten Clustergrenzen und Fortsetzungswerte bzw. später die Regressionskoeffizienten. Definiere die Zuordnung neuer Zustände einschließlich Randfällen außerhalb der Trainingsintervalle. Wende diese fixierten Regeln auf den Bewertungspfaden in chronologischer Reihenfolge an. Vergleiche die tatsächlich realisierten Ausübungen getrennt von hypothetischen Entscheidungen in der Rückwärtsinduktion.

Für gemeinsame Cluster ist die Zuordnung kontrolliert zu halten; für getrennte Cluster darf jede Variante ihre eigene trainierte Policy anwenden. Dokumentiere den Unterschied. Der Out-of-Sample-Wert einer zulässigen Policy ist im Erwartungswert eine Untergrenze für den optimalen Wert im jeweils verwendeten Simulationsmodell, nicht notwendigerweise jeder einzelne Monte-Carlo-Schätzwert für das kontinuierliche Modell.

## 13. Robustheits- und Sensitivitätsuntersuchungen

Die Anwendung soll so strukturiert sein, dass der Vergleich unkompliziert für verschiedene numerische Einstellungen wiederholt werden kann.

Mindestens folgende Größen sollen konfigurierbar sein:

- Anzahl der Monte-Carlo-Pfade,
- Anzahl beziehungsweise Größe der Zeitschritte,
- Anzahl der Quantilcluster,
- Zufallsseed,
- `MathContext` der `BigDecimal`-Variante,
- später gegebenenfalls die Longstaff-Schwartz-Regressionsbasis.

Damit sollen insbesondere folgende Experimente möglich sein:

1. steigende Anzahl von Monte-Carlo-Pfaden,
2. feinere Euler-Zeitdiskretisierung,
3. unterschiedliche Anzahl von Quantilclustern,
4. unterschiedliche `BigDecimal`-Präzisionen,
5. gemeinsame gegenüber getrennten Quantilclustern,
6. Quantilcluster gegenüber Longstaff-Schwartz.

Die Architektur soll diese Vergleiche ermöglichen, ohne dafür ein allgemeines Simulationsframework zu entwickeln.

Ergänze im Laufe des Projekts gezielte Parameter- und Lastszenarien:

- unterschiedliche Laufzeiten, Strikes (ITM/ATM/OTM), Nominale und Anzahl Exercise Dates,
- niedrige/hohe Zinsen, Volatilität und Mean Reversion,
- Werte nahe der Feller-Grenze und später explizit zugelassene Feller-Verletzungen,
- Präzisionen 16, 20, 34, 50, 64 und punktuell 100 Stellen sowie unterschiedliche Rundungsmodi,
- gezielt konstruierte Fälle nahe einer Clustergrenze oder Exercise Boundary,
- unterschiedliche Portfolio- und Threadgrößen.

Variiere jeweils möglichst nur einen Faktor und vermeide das vollständige kartesische Produkt. Unterscheide einen schnellen Smoke-Test, ein Standardexperiment und explizit gestartete umfangreiche Research-Läufe.

## 14. Tests

Es sollen automatisierte Tests für fachlich und numerisch relevante Eigenschaften vorgesehen werden.

Dazu gehören mindestens:

- Ablehnung von CIR-Parametern, die die Feller-Bedingung verletzen,
- Reproduzierbarkeit bei festem Seed,
- identische Zufallszahlen für beide numerischen Varianten,
- korrekte CIR-Zerobondpreise für ausgewählte Testfälle,
- korrekte Cashflows der Kuponanleihe,
- korrekte Immediate Exercise Values,
- sinnvolles Verhalten der Diskontierung,
- korrekte Bildung gleich großer beziehungsweise nahezu gleich großer Quantilcluster,
- deterministisches Verhalten bei Gleichständen der Clustervariable,
- korrekte Berechnung der Cluster-Mittelwerte,
- Konsistenz der Rückwärtsinduktion,
- nachvollziehbare Exercise Decisions in bewusst konstruierten einfachen Fällen,
- korrekte Aktualisierung des zukünftigen optimalen Cashflows nach einer früheren Exercise Decision,
- Vergleich von `double` und `BigDecimal` in einfachen Situationen, in denen nahezu identische Ergebnisse erwartet werden.

Tests sollen deterministisch sein. Monte-Carlo-Tests sollen nicht aufgrund zufälliger Schwankungen sporadisch fehlschlagen.

Ergänze die Tests im jeweiligen Entwicklungsschritt durch:

- unabhängige Referenzwerte für CIR-Zerobondpreise und später europäische Optionspreise; keine Referenzwerte ausschließlich aus derselben Implementierung erzeugen,
- Pfad-Cashflows ohne Ausübung, Exercise-Entscheidung bei Gleichstand, Cashflows genau am Exercise Date,
- gemeinsamer Cluster-Modus mit identischen Pfad-IDs und getrennte Cluster mit gezielt konstruiertem Wechsel,
- korrekte Zuordnung unabhängiger Bewertungspfade zu trainierten Clustern,
- Reduktion der bermudischen auf die europäische Option bei nur einem Exercise Date und identischem Underlying,
- Nichtnegativität, monotones Verhalten bezüglich Strike und lineare Skalierung bei gemeinsamer Skalierung von Cashflows und Strike,
- Plausibilitätsprüfung, dass der theoretische bermudische Wert nicht unter dem einer ansonsten identischen europäischen Option an einem zulässigen Exercise Date liegt (Monte-Carlo-Schätzungen und approximative Policies können davon abweichen),
- Stabilität der höherpräzisen Referenz und nachvollziehbare statistische Toleranzen.

## 15. Technische Rahmenbedingungen

Das Projekt soll als Maven-Projekt in Java 21 umgesetzt werden. Projekt- und Artefaktname: `cir-numerical-precision-lab`.

Verwende JUnit 5, UTF-8 und Git. Wähle tatsächlich verfügbare, kompatible Versionen der Abhängigkeiten. Verwende kein Spring Boot, Quarkus, keine Datenbank und kein allgemeines Simulationsframework. Quellcode, Bezeichner, Kommentare und technische Projektdokumentation sollen auf Englisch verfasst werden.

Java 8 ist die bisherige Laufzeitumgebung unserer mathematischen Modelle. Ein späterer Java-8/21/25-Vergleich ist daher interessant, aber keine Voraussetzung für den Java-21-Prototyp. Falls er durchgeführt wird, verwende für den isolierten Laufzeitvergleich einen gemeinsamen Java-8-kompatiblen Rechenkern oder ein separates kompatibles Benchmark-Modul mit denselben Algorithmen und Bibliotheken. Vergleiche gegebenenfalls `Math` und `StrictMath` und dokumentiere Floating-Point-Semantik, `strictfp`, JVM-Version, JIT und Garbage Collector. Eine Java-21-Anwendung selbst ist nicht auf Java 8 ausführbar.

Plane zunächst eine einfache Maven-Struktur mit Komponenten für Modell und Cashflows, analytische Bondbewertung, Zufallszahlen, Simulation, Discounting, numerische Varianten, Continuation Estimation, Backward Induction, Experimente, Statistik und Reporting. Halte gemeinsame fachliche Konventionen und Testfälle zentral, aber implementiere die rechenintensiven Double- und BigDecimal-Algorithmen separat und nachvollziehbar.

Vermeide unnötige Frameworks, Abstraktionen und Dependencies.

Es handelt sich ausdrücklich um einen kleinen mathematisch-numerischen Prototypen und nicht um eine produktive Anwendung.

Die Architektur soll so einfach wie möglich bleiben, aber dennoch sauber genug sein, um

- `double` und `BigDecimal`,
- gemeinsame und getrennte Quantilcluster sowie
- später Quantilclustering und Longstaff-Schwartz

systematisch miteinander vergleichen zu können.

Insbesondere soll vermieden werden, durch übermäßige Abstraktion beide numerischen Implementierungen künstlich in eine gemeinsame generische Numeric-API zu pressen, wenn dadurch die Implementierung unnötig kompliziert oder der eigentliche Vergleich weniger transparent wird.

## 16. Erwartete Fragestellungen des Experiments

Das Projekt soll insbesondere folgende Hypothesen empirisch untersuchen:

1. Die direkten arithmetischen Abweichungen zwischen `double` und `BigDecimal` sind bei typischen Modellparametern sehr klein.

2. Kleine Abweichungen können dennoch nahe einer Clustergrenze oder Exercise Boundary zu diskret unterschiedlichen Entscheidungen führen.

3. Durch die Rückwärtsinduktion können sich solche Unterschiede auf frühere Exercise Dates fortpflanzen.

4. Trotz dieser Effekte bleibt der Unterschied zwischen `double` und `BigDecimal` möglicherweise klein gegenüber
   - Monte-Carlo-Fehler,
   - Euler-Diskretisierungsfehler und
   - Approximationsfehler bei der Schätzung des Continuation Value.

5. Die Wahl der Methode zur Schätzung des Continuation Value, insbesondere Quantilclustering gegenüber Longstaff-Schwartz, könnte einen größeren Einfluss auf den Optionspreis und die Exercise Policy besitzen als die Wahl zwischen `double` und `BigDecimal`.

6. Der Performance- und Speichermehraufwand höherer Präzision könnte insbesondere bei vielen Pfaden, Exercise Dates und Verträgen relevant sein.

7. Kompensierte Double-Summation oder gezielte hybride Verfahren könnten einen Teil eines beobachteten Präzisionsvorteils zu geringeren Kosten erzielen.

8. Die Wahl des Diskretisierungsverfahrens und die Güte der geschätzten Exercise Policy könnten wichtiger sein als eine weitere Erhöhung der Arbeitspräzision. Unabhängige Referenzen sollen dies untersuchbar machen.

Diese Hypothesen sollen ausdrücklich nicht als vorgegebenes Ergebnis behandelt werden. Die Anwendung soll so gestaltet werden, dass sie möglichst transparent empirisch überprüft werden können.


## 17. Analytische Referenzbewertung und fachliche Validierung

Die folgenden Referenzfälle ergänzen das Hauptprojekt. Sie sollen die Richtigkeit der gemeinsamen Modellkomponenten und des Monte-Carlo-Verfahrens prüfen, nicht die bermudische Option als Forschungsgegenstand ersetzen.

### 17.1 Analytischer CIR-Zerobond

Implementiere und validiere zunächst in beiden Arithmetiken die analytische CIR-Zerobondformel

$$
P(t,T)=A(t,T)\exp(-B(t,T)r_t).
$$

Sie wird bereits für die Immediate Exercise Values der bermudischen Option benötigt. Dokumentiere die Formeln in der in Abschnitt 1 verwendeten Parameternotation und validiere sie gegen unabhängige Testwerte. Prüfe numerisch sensible Fälle, insbesondere kurze Restlaufzeiten, kleine Volatilitäten und große Exponenten. Berechne die Bondpreise nicht zusätzlich über Monte Carlo.

### 17.2 Europäische Option auf einen Zerobond

Verwende als Validierungsexperiment einen europäischen Call (optional zusätzlich Put) auf einen Zerobond mit Optionsfälligkeit $T$ und Bondfälligkeit $S>T$. Für einen Einheitsbond und einen Strike $K$ gilt:

$$
V_0=\mathbb E^Q\left[\exp\left(-\int_0^T r_s ds\right)(P(T,S)-K)^+\right].
$$

Skaliere bei einem nicht normierten Vertrag die Cashflows und den Strike konsistent. Vergleiche den Monte-Carlo-Schätzer mit der analytischen CIR-Optionsformel unter Verwendung der nichtzentralen Chi-Quadrat-Verteilung. Die Simulation soll denselben CIR-Pfadgenerator und dieselbe Diskontierung wie die bermudische Bewertung nutzen. Prüfe optional die Put-Call-Parität $C_0-P_0=P(0,S)-KP(0,T)$ für den Einheitsbond.

Eine vollständig hochpräzise Implementierung der nichtzentralen Chi-Quadrat-CDF kann anspruchsvoll sein. Prüfe geeignete Bibliotheken bzw. etablierte Algorithmen (z. B. stabil ausgewertete Poisson-Mischung regulierter Gammafunktionen). Berücksichtige Konvergenz, Abschneidefehler, Verteilungsschwänze, Underflow/Overflow und große Nichtzentralitätsparameter. Validiere gegen unabhängige Referenzwerte aus R, Python oder anderer geeigneter Software und dokumentiere den validierten Parameterbereich.

Priorisiere zunächst eine unabhängig validierte analytische Referenz in geeigneter Genauigkeit. Die vollständige BigDecimal-CDF ist ein eigener späterer Arbeitsschritt, falls sie den Prototyp verzögern würde. Kennzeichne eine Double-basierte Referenz eindeutig; gib sie nicht als hochpräzise BigDecimal-Lösung aus. Keine stillschweigende Ersatzimplementierung und keine erfundenen Referenzwerte.

### 17.3 Europäische Option auf eine Kuponanleihe

Eine europäische Option auf den allgemeinen positiven Cashflow-Strom ist der natürliche Kontrollfall für die bermudische Option mit genau einem Exercise Date. Prüfe später die Jamshidian-Zerlegung im Einfaktor-CIR-Modell. Dokumentiere deren Voraussetzungen (insbesondere Monotonie und gegebenenfalls Existenz eines kritischen Zinssatzes), bestimme die zugehörigen Einzelstrikes mit einer geeigneten Nullstellensuche und behandle Randfälle ohne innere Nullstelle explizit.

Vergleiche die Referenzbewertung mit beiden Monte-Carlo-Varianten und dem Ein-Exercise-Date-Spezialfall der bermudischen Implementierung. Eine analytische Formel für den allgemeinen bermudischen Vertrag wird nicht vorausgesetzt.

### 17.4 Optionale unabhängige bermudische Referenz

Als spätere Erweiterung prüfe eine deterministische Bewertung mittels Finite-Differenzen-Verfahren oder eines geeigneten eindimensionalen Gitterverfahrens. Die CIR-Bewertungsgleichung zwischen den Exercise Dates lautet in unserer Notation:

$$
\partial_t V+\theta(\mu-r)\partial_r V+\frac{1}{2}\sigma^2 r\partial_{rr}V-rV=0.
$$

An den Exercise Dates ist die Maximumsbedingung mit dem Immediate Exercise Value anzuwenden. Randbedingungen, Zinsgitter, Zeitdiskretisierung und Gitterkonvergenz sind fachlich zu begründen. Diese Referenz ist unabhängig von der Monte-Carlo-Fortsetzungswertschätzung, aber selbst nicht frei von numerischen Fehlern. Sie ist kein Bestandteil des ersten Prototyps.

## 18. Performance und Portfoliobewertung

### 18.1 Laufzeit und Ressourcen

Zeige schon im ersten Prototyp die Berechnungsdauer je Rechenvariante und den Preisvergleich. Eine einfache Messung mit `System.nanoTime()` ist hierfür als orientierende Zeitmessung ausreichend, jedoch nicht als alleinige Grundlage belastbarer Performance-Aussagen.

Ergänze später JMH für reproduzierbare Benchmarks. Berücksichtige Warm-up, mehrere Messiterationen und Forks, Dead Code Elimination, Constant Folding, Objektallokationen und Garbage Collection. Trenne die folgenden Messbereiche, soweit sinnvoll:

- analytische Bond- und Optionsbewertung,
- Pfadsimulation und Diskontierung,
- Berechnung der Immediate Exercise Values,
- Sortierung und Quantilclusterbildung,
- Mittelwerte bzw. später Regression und Backward Induction,
- vollständige Bewertung (End-to-End) und Aggregation.

Bei isolierten Bewertungsbenchmarks dürfen die gemeinsamen Zufallszahlen außerhalb der Messung vorbereitet werden; bei End-to-End-Messungen sind deren Erzeugung und erforderliche Konvertierungen einzubeziehen. Ordne gemeinsame Vorberechnungen transparent zu. Logging, CSV und HTML-Erstellung gehören nicht zur gemessenen Bewertungszeit. Messe auch ohne aktivierte Detaildiagnostik. Vergleiche vergleichbare fachliche Arbeit und dokumentiere Abweichungen.

Erfasse Zeit, Throughput, Speedup/Overhead-Faktor, allokierte Bytes je Operation, GC-Statistiken sowie bei Bedarf Peak Heap Memory. Unterscheide Speicherallokation von maximalem Speicherbedarf. Gib auch die verwendete Java-Version, JVM-Argumente, CPU, Betriebssystem, Heap-Konfiguration und Bibliotheksversionen an. Halte die Messbedingungen zwischen den Varianten möglichst konstant.

### 18.2 Einzelvertrag und Vertragsbestand

Unterstütze schrittweise Portfolios bermudischer Optionen mit unterschiedlichen Cashflows, Strikes, Nominalen und Exercise Dates. Erlaube auch europäische Referenzverträge. Unterscheide homogene Lasttests und heterogene Portfolios.

Untersuche zwei Strategien separat:

1. Independent valuation: eigenständige Bewertung jedes Vertrags.
2. Shared market simulation: Wiederverwendung kompatibler CIR-Pfade und ggf. fachlich zulässiger Zwischenergebnisse für mehrere Verträge bei identischem Modell und passendem Zeitraster. Die Exercise Policies und Cashflows bleiben vertragsabhängig.

Erfasse Gesamtzeit, durchschnittliche Zeit je Vertrag, Verträge pro Sekunde, Speicherbedarf und Differenz des aggregierten Portfoliowerts. Vermeide versehentliche Cache-Treffer, die nur bereits berechnete Vertragswerte zurückliefern. Bei gemeinsamen Pfaden sind die Vertragswerte statistisch abhängig; bestimme einen Portfolio-Standardfehler gegebenenfalls aus den pfadweise aggregierten Cashflows, nicht durch Addition einzelner Standardfehler unter unbegründeter Unabhängigkeitsannahme.

Vorgeschlagene spätere Lastreihen (nicht als vollständiges kartesisches Produkt):

| Faktor | Beispielwerte |
| --- | --- |
| Verträge | 1, 10, 100, 1.000, 10.000 |
| Trainings-/Bewertungspfade je Menge | 100, 1.000, 10.000, 100.000 |
| Zeitschritte pro Jahr | 12, 52, 252 |
| Quantilcluster | 5, 10, 20, 50 |
| BigDecimal-Präzision | 16, 20, 34, 64 |
| Threads | 1, 2, 4, 8 |

Beginne sequenziell. Parallelisierung ist eine spätere Erweiterung; verwende kontrollierbare Threadzahlen, deterministische Pfadzuordnung und eine definierte Aggregationsreihenfolge. Sehr große BigDecimal-Experimente sollen bewusst gestartet werden, nicht Bestandteil des Standardlaufs sein.

## 19. Konfiguration, Ergebnisdarstellung und Dokumentation

### 19.1 Erster Testfall

Schlage im Entwicklungsplan einen kleinen, vollständig spezifizierten und fachlich plausiblen Testfall vor. Als Ausgangspunkt können diese Werte dienen (keine bereits validierten Preisreferenzen):

```properties
cir.initialRate=0.03
cir.meanReversionSpeed=0.50
cir.longTermMean=0.04
cir.volatility=0.10
bond.notional=100.00
bond.couponRate=0.04
bond.maturity=5.0
bond.couponFrequency=1
option.strike=100.00
option.exerciseDates=1.0,2.0,3.0,4.0
simulation.stepsPerYear=12
simulation.paths=1000
simulation.seed=42
cluster.count=10
decimal.precision=34
decimal.roundingMode=HALF_EVEN
```

Die Modellparameter erfüllen die Feller-Bedingung. Die Zeitkonvention für Cashflows an Exercise Dates aus Abschnitt 2 gilt auch hier. Prüfe, ob dieser Strike im Zusammenspiel mit den verbleibenden Cashflows genügend ITM-Pfade und unterschiedliche Ausübungen erzeugt; ergänze bei Bedarf einen geeigneteren Testfall. Für Out-of-Sample-Experimente sind separate Pfadanzahlen und Seeds für Training und Bewertung zu konfigurieren.

### 19.2 Ergebnisse

Der erste Prototyp soll eine verständliche Konsolenausgabe liefern: Preis je Variante, absolute und relative Differenz, Standardfehler, Berechnungsdauer sowie zentrale Cluster- und Exercise-Metriken. Gib den Namen der Methode, den Clustermodus und die BigDecimal-Präzision an.

Ergänze CSV-Exporte für die wissenschaftliche Auswertung in R, Python oder Excel. Trenne bei Bedarf Dateien für Experimentzusammenfassungen, Exercise-Date-/Cluster-Statistiken, ausgewählte Pfadverläufe und Benchmarks. Schreibe hochpräzise Werte als Dezimalstrings ohne vorherige Double-Konvertierung.

Speichere pro Lauf eine eindeutige Experiment-ID, vollständige Konfiguration, Modell- und Vertragsparameter, Seeds, Pfadanzahlen, Zeitraster, Verfahren, Clustermodus, Präzision, Rundungsmodus, Laufzeitumgebung und Git-Commit (soweit verfügbar). Unterscheide Training und Bewertung, Diagnose- und Benchmark-Läufe sowie analytische Referenz und höherpräzise Algorithmusreferenz.

Ein späterer offline nutzbarer HTML-Report kann folgende Diagramme erzeugen:

- Optionspreis und Konfidenzintervall in Abhängigkeit von Pfadzahl und Zeitschrittgröße,
- numerische Differenz in Abhängigkeit von Präzision und Simulationszeit,
- Clusterwechsel, Unterschiede der Continuation Values und Exercise Rates je Exercise Date,
- Abstand zur Exercise Boundary und Unterschiede zwischen den Methoden,
- Berechnungsdauer und Throughput in Abhängigkeit von Portfolio- und Pfadgröße,
- Speicherallokation und optional Thread-Skalierung.

Verwende bei Bedarf logarithmische Achsen und kennzeichne Nullwerte angemessen. Keine erfundenen Messwerte, keine unbegründeten Schlussfolgerungen. Ein Dashboard ist kein Muss für den ersten Prototyp.

### 19.3 Wissenschaftliche Dokumentation

Erstelle im späteren Projekt eine README mit Build- und Ausführungsanweisungen, eine mathematische Beschreibung der Formeln und Quellen sowie eine Dokumentation der Methodik, Experimente und bekannten Grenzen. Verifiziere wissenschaftliche Referenzen, Bibliotheken und externe Vergleichswerte. Dokumentiere insbesondere unterschiedliche Fehlerquellen, Rundungsentscheidungen, Zufallszahlengrenzen, Pfadauswahl, Feller- und Euler-Konvention, Clusterzuordnung, Trainings-/Bewertungsprotokoll und Messmethodik.

## 20. Vorgehen

Ich möchte das Projekt ausdrücklich schrittweise mit dir entwickeln.

**Erstelle deshalb jetzt noch keinen Code und keine Dateien und führe keinerlei Änderungen am Dateisystem durch.**

Analysiere zunächst ausschließlich die Anforderungen und erstelle einen sinnvollen Entwicklungsplan.

Der Plan soll insbesondere behandeln:

1. fachliche Modellierung von CIR-Modell, Bond und bermudischer Option,
2. mathematische Formulierung der Rückwärtsinduktion,
3. genaue Definition der Quantilcluster-Methode,
4. Entscheidung, welche Pfade in die Cluster-Schätzung eingehen,
5. Strategie für gemeinsame und getrennte Cluster,
6. vorgesehene Java-Komponenten und Verantwortlichkeiten,
7. Zusammenspiel von `double`- und `BigDecimal`-Implementierung,
8. Strategie für Zufallszahlen und Common Random Numbers,
9. Zeitdiskretisierung und Diskontierung,
10. Umgang mit möglichen negativen Euler-Werten,
11. Definition und Erhebung der beschriebenen Vergleichsmetriken,
12. spätere Ergänzung von Longstaff-Schwartz,
13. vorgesehene automatisierte Tests,
14. sinnvolle Standardparameter für einen ersten kleinen Testlauf,
15. schrittweise Reihenfolge der späteren Implementierung,
16. Strategie für unabhängige Validierung mittels europäischer Referenzfälle und ggf. einer späteren bermudischen Referenz,
17. Trennung von Trainings- und Bewertungspfaden sowie gepaarte statistische Analyse,
18. geplante Präzisionsstrategie, geeignete Bibliotheken und deren Grenzen,
19. Laufzeitmessung, JMH, Portfoliobewertung und spätere Skalierung,
20. Auswertung, Dokumentation und klare Abgrenzung zwischen MVP und optionalen Erweiterungen.

Schlage eine inkrementelle Entwicklung in etwa folgender Reihenfolge vor:

| Phase | Inhalt | Priorität |
| --- | --- | --- |
| 0 | Anforderungen, Mathematik, Architektur, Validierungsstrategie | Jetzt: ausschließlich Planung |
| 1 | CIR-Modell, Cashflow-Modell, Double-/BigDecimal-Zerobondpreise, gemeinsame Zufallszahlen | Grundlage |
| 2 | CIR-Simulation, Diskontierung, europäischer Ein-Exercise-Date-MC-Kontrollfall und unabhängige Preisreferenz | Grundlage / Validierung |
| 3 | Bermudische Rückwärtsinduktion mit Quantilclustern, gemeinsame/getrennte Cluster, Detailmetriken, Konsolenausgabe und erste Laufzeitmessung | Kern des MVP |
| 4 | Out-of-Sample-Policies, gepaarte Statistik, Sensitivität und Präzisionsreihen, CSV | Belastbare Experimente |
| 5 | Longstaff-Schwartz mit geeigneter stabiler Regression | Fachliche Erweiterung |
| 6 | JMH und Portfolios, Speicheranalyse, ggf. Parallelisierung | Performance und Skalierung |
| 7 | Jamshidian-Referenz, höherpräzise analytische CDF soweit noch offen, weitere Diskretisierung, deterministische bermudische Referenz, Hybrid- und JVM-Vergleiche, HTML | Optionale Vertiefung |

Die Phasen dürfen fachlich begründet umgeordnet werden. Die schwierige hochpräzise analytische Optionsformel darf insbesondere nicht den bermudischen Kernprototyp blockieren. Eine verifizierte analytische Referenz für den europäischen Kontrollfall ist dennoch ein wichtiges Validierungsziel. Nach späterer Freigabe jeder Implementierungsphase: Build und Tests ausführen, Entscheidungen und offene Punkte dokumentieren, keine Platzhalter als fertig deklarieren und keine Messergebnisse erfinden. Kein automatischer Push in ein externes Repository.

Erkläre in deinem Plan auch, welche fachlichen Konventionen bereits festgelegt sind und welche Entscheidungen noch zur gemeinsamen Abstimmung stehen. Nenne konkrete Akzeptanzkriterien für den MVP, sinnvolle Standardparameter und Risiken (insbesondere Speicherbedarf, Schätzung der Fortsetzungswerte, Vergleichbarkeit der Arithmetik und Stabilität bei diskreten Entscheidungen).

Weise im Plan ausdrücklich auf Stellen hin, an denen der gewünschte `double`-gegen-`BigDecimal`-Vergleich methodisch problematisch oder möglicherweise nicht aussagekräftig ist.

Vermeide unnötige Komplexität. Falls Anforderungen noch unklar sind, triff für einen Prototypen sinnvolle Annahmen und kennzeichne diese ausdrücklich.

**Wichtig: Führe noch keine Änderungen durch. Gib mir in deiner ersten Antwort ausschließlich deine Analyse und deinen Entwicklungsplan.**
