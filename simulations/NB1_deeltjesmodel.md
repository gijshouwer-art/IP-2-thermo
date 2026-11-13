# NB 1: Deeltjesmodel


Computer simulaties zijn een belangrijk onderdeel bij natuurkundig onderzoek: wat we niet in het lab kunnen doen, of waar we niet goed aan kunnen meten, kunnen we alsnog bestuderen dankzij simulaties. In deze module leggen we de basis voor natuurkundige computer simulaties en maken we de eerste stappen in ons deeltjes model.

De basis van een goede computer simulatie ligt, net als bij het oplossen van een natuurkundig probleem in een gedegen systematische aanpak. Deze begint vaak met:
1. Het maken van een tekening van de situatie
2. Het noteren van de werkformules
3. Het aangeven van richtingen en noteren van variabelen. 

```{exercise}
Als startopdracht modelleren we twee fietsers. Doorloop de bovenstaande drie stappen. 
```

<!-- Jouw antwoord hier -->
**1. Tekening.** Een rechte weg (de $x$-as) met de oorsprong bij de start. Beide fietsers starten op $t=0$ in $x=0$ en rijden in de positieve $x$-richting (naar rechts):
- fietser 1 rijdt met een constante snelheid $v_1$ (pijl van vaste lengte);
- fietser 2 start vanuit stilstand en heeft een constante versnelling $a_2$ (pijlen die steeds langer worden).

**2. Werkformules.**
- Fietser 1 (eenparige beweging): $x_1(t) = x_{1,0} + v_1 t$
- Fietser 2 (eenparig versnelde beweging): $v_2(t) = v_{2,0} + a_2 t$ en $x_2(t) = x_{2,0} + v_{2,0}t + \tfrac{1}{2}a_2t^2$
- Numeriek (per tijdstap $\Delta t$): $x_{i+1} = x_i + v_i\Delta t$ en $v_{i+1} = v_i + a\Delta t$

**3. Richtingen en variabelen.** Positief is naar rechts.
| grootheid | symbool | waarde |
|---|---|---|
| tijdstap | $\Delta t$ | 0,25 s |
| eindtijd | $t_e$ | 10 s |
| beginposities | $x_{1,0}, x_{2,0}$ | 0 m |
| snelheid fietser 1 | $v_1$ | 5 m/s |
| beginsnelheid fietser 2 | $v_{2,0}$ | 0 m/s |
| versnelling fietser 2 | $a_2$ | 0,5 m/s² |




## Numerieke modellen
Om een simulatie te maken, splitsen we het volledige proces op in kleine stappen met een vaste (tijd)interval, zoals in onderstaande figuur. Elk stapje krijgt een nummer, van $0$ tot $N$. We gebruiken de letter $i$ om naar een specifieke stap te verwijzen.

De gegevens in stap $i$ worden gebruikt om de waarden in stap $i+1$ te berekenen, bijvoorbeeld:

$$x_{i+1} = x_i + v_i\cdot\Delta t$$

```{figure} ../figures/interval.jpg
:width: 50%

Het gehele proces wordt opgebroken in kleine stukken: tijdsdiscretisatie.
```
Wellicht heb je bovenstaande al eens gezien met het programma Coach. Anders dan bij Coach moeten we expliciet de waarden opslaan in een array. Omdat het steeds aanmaken van een nieuwe cel in een array veel rekentijd kost, is het handiger om eerst arrays te maken met waarde 0.

We nemen als voorbeeld twee fietsers, waarbij de fietser 1 een eenparige beweging uitvoert en fietser 2 een eenparig versnelde beweging:

```python
import numpy as np
import matplotlib.pyplot as plt

# Beginwaarden
dt = 0.25                  #s
t_e = 10                   #s
N = int(t_e/dt)+1
t = np.arange(0,t_e+dt,dt) #s

# fietser 1
x_1 = np.zeros(N)          #m     maakt array van lengte N met alleen maar waarde 0.
v_1 = 5                    #m/s

# fietser 2
x_2 = np.zeros(N)          #m     
v_2 = np.zeros(N)          #m/s     
a_2 = .5                   #m/s^2     

# Uitvoeren van de simulatie
for i in range(0,N-1):
    x_1[i+1] = x_1[i] + v_1*dt
    
    x_2[i+1] = x_2[i] + v_2[i]*dt
    v_2[i+1] = v_2[i] + a_2*dt

# Plotten van de gegevens.    
plt.figure()

plt.xlabel('$t$(s)')
plt.ylabel('$x$(m)')

plt.plot(t,x_1,'r.', markersize=1, label='fietser 1')
plt.plot(t,x_2,'b.', markersize=1, label='fietser 2')

plt.show()
print(x_2)
print(x_1)
```

Wat op zou kunnen vallen is dat de snelheid in een cel constant gehouden wordt ($x_i = v_i\cdot\Delta t$), terwijl de snelheid tijdens het tijdsinterval $\Delta t$ toeneemt. De oplossing is dus een benadering van de exacte oplossing, maar door deze methode (Euler explicit forward), moeten we wel voorzichtig zijn met de grootte van het tijdsinterval.

```{exercise}
Maak het tijdsinterval $dt$ groter en vergelijk de afgelegde weg in de tien seconde. Wat gebeurt er met de afgelegde weg als we $dt$ groter maken? Is dat logisch vanuit gezien de natuurkunde? En gezien de methode van discretisatie? Leg uit.
```


<!-- Jouw antwoord hier  -->
Afgelegde weg van fietser 2 na 10 s (exact: $\tfrac{1}{2}a_2t^2 = 25$ m):

| $\Delta t$ (s) | 0,01 | 0,25 | 0,5 | 1 | 2 |
|---|---|---|---|---|---|
| $x_2(10)$ (m) | 24,975 | 24,375 | 23,75 | 22,5 | 20,0 |

Voor fietser 1 verandert er niets: die blijft op 50 m uitkomen, want zijn snelheid is constant en dus maakt de grootte van de stap niet uit.

Voor fietser 2 wordt de afgelegde weg *kleiner* als $\Delta t$ groter wordt. Natuurkundig is dat niet logisch: de beweging zelf hangt niet af van hoe wij de tijd opdelen, de echte afstand is altijd 25 m. Vanuit de discretisatie is het wel logisch: in elke stap gebruiken we de snelheid $v_i$ aan het *begin* van het interval ($x_{i+1} = x_i + v_i\Delta t$), terwijl de snelheid tijdens dat interval toeneemt. We onderschatten dus elke stap de verplaatsing. Opgeteld geeft dat $x_2 = \tfrac{1}{2}a_2t_e^2 - \tfrac{1}{2}a_2t_e\Delta t$: de fout is evenredig met $\Delta t$ (Euler forward is een eerste-orde methode). Een kleinere $\Delta t$ geeft een nauwkeuriger antwoord, maar kost meer rekenstappen.

```python
  #clears memory as we don't need it further down  
%reset -f 
```

## Naar een deeltjesmodel
In Q2 werken we aan een deeltjes model. Dat model bestaat uiteindelijk uit een stuk of 100 deeltjes die in een afgesloten volume bewegen. Om voor elk deeltje apart de eigenschappen te bepalen, wordt enorm veel code. Omdat de deeltjes dezelfde eigenschappen hebben (massa, radius, snelheid, positie) kunnen we gebruik maken van classes.

Een *class* is een blauwdruk voor het creëren van objecten (zoals deeltjes, atomen, planeten, auto's, studenten... noem maar op!). Het stelt je in staat om data (zoals positie, massa) en gedrag (zoals bewegen, botsen) te bundelen tot één overzichtelijk geheel. Precies wat we nodig hebben. Laten we de anatomie van een deeltje creëren.

```python
# Definieer een class voor een deeltje. 
# Een class beschrijft alleen het type eigenschappen dat een ding heeft, bijvoorbeeld: in de beschrijving van de class zeg je:
# - een deeltje heeft een positie.
# - een deeltje heeft een functie die, wanneer aangeroepen, de positie van dat deeltje bijwerkt.
# Merk op dat de class zelf geen deeltje is! 
# Een object dat tot een bepaalde class behoort, wordt een instant van die class genoemd. 
# Je kunt meerdere instants van een class hebben.
import numpy as np
import matplotlib.pyplot as plt

class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m                         # massa van het deeltje
        self.v = np.array(v, dtype=float)  # snelheids vector
        self.r = np.array(r, dtype=float)  # positie vector
        self.R = np.array(R, dtype=float)  # radius van het deeltje

    def update_position(self,dt):
        """Werk de positie van het deeltje bij op basis van zijn snelheid en tijdstap dt."""
        self.r += self.v * dt

print(ParticleClass)
```

We hebben alleen de class aangemaakt, maar nog niet een object zelf. Dat moeten we dus nog doen. 

We hebben ook met onze class de mogelijkheid gemaakt om een nieuwe positie uit te rekenen: de volgende positie is de huidige positie + de snelheid maal de tijdstap, ofwel: 

$$ x_{n+1}=x_n + v\cdot\Delta t $$

Laten we nu een ook echt een deeltje creëren en kijken wat we ermee kunnen doen.

```python
# een enkel deeltje:
this_particle = ParticleClass(m=1.0, v=[5.0, 0], r=[0.0, 0.0], R=1.0)
print(this_particle)
```

```python
%whos
```

Nu we een object hebben, laten we eens kijken of we de eigenschappen van het object kunnen opvragen:

```python
print("Massa van het deeltje " + str(this_particle.m) + " kg")
print("Huidige positie is " + str(this_particle.r))

# oproepen van de functie 'update_position' met dt = 1.0 seconde
this_particle.update_position(1.0)
print("Nieuwe positie " + str(this_particle.r))
```

Let op! Run je bovenstaande cell opnieuw, dan zal de positie ook steeds veranderen!

Laten we nog een deeltje maken om te laten zien dat het aanroepen van update op het ene deeltje geen invloed heeft op het andere deeltje: het zijn onafhankelijke objecten!

```python
that_particle = ParticleClass(m=1.0, v=[2.0, 0], r=[5.0, 0.0],R=1.0)
```

```python
print("Huidige positie van this_particle " + str(this_particle.r))
print("Huidige position van that_particle " + str(that_particle.r))

# oproepen van de functie 'update_position' met dt = 2.0 seconde
that_particle.update_position(2.0)
print("Nieuwe positie van this_particle " + str(this_particle.r))
print("Nieuwe positie van that_particle " + str(that_particle.r))
```

We kunnen ook een array van objecten maken.
Dat zou eenvoudig moeten kunnen als het gaat om dezelfde deeltjes (zelfde massa, snelheid, maar andere startpositie). 

```python
particle_array = []

for start_x in range(4):
    particle_array.append(ParticleClass(m=1.0, v=[2.0, 0], r=[start_x, 0.0],R=1.0))

print("Hoe het array eruit ziet volgens python: " + str(particle_array))
print("Hoe een element van het array eruit ziet volgens python: " + str(particle_array[0]))
print("Hoe een eigenschap (in dit geval: 'r') van het element van het array eruit ziet volgens python: " + str(particle_array[0].r))
```

```{exercise} Voorspel
:label: ex-intro-1

Gegeven bovenstaande particles die je hebt aangemaakt, hoe zien de posities er dan uit van die deeltjes?
Voorspel en controleer vervolgens met onderstaande code.
```

```python
import numpy as np
import matplotlib.pyplot as plt

plt.figure()
for particle in particle_array:
    plt.plot(particle.r[0],particle.r[1],'k.')
plt.show()

```

Wat eerder op moet zijn gevallen, en wat hierboven ook weer duidelijk wordt is dat de positie een 2D vector is.
We kunnen dus onafhankelijk in de $x$ en $y$ richting bewegen!

Door gebruik te maken van een array waarin alle deeltjes zijn opgeslagen, kun je voor elk deeltje dezelfde bewerking uitvoeren, bijvoorbeeld allemaal een stukje updaten in de tijd!
Hierbij kan je mooi gebruik maken van hoe Python loops maakt: als je één voor één de verschillende elementen van array af wilt gaan en daar dezelfde bewerking op doen, gaat dat zo:

```python
for particle, particle_object in enumerate(particle_array):
    print("Deeltje " + str(particle) + ", huidige positie: " + str(particle_object.r))
    print("oproepen van de functie 'update_position' met dt = 1.0 second ")
    particle_object.update_position(1.0)
    print("Volgende positie " + str(particle_object.r))
    if particle < len(particle_array) - 1:
        print("Naar het volgende deeltje \n")
    
```

Let op! In bovenstaande code maken we handig gebruik van twee programmeerconcepten.
`enumerate` hangt een nummer (counter) aan elk item in de array. Zo kun je bijhouden met welk deeltje je bezig bent. 

Omdat onze code steeds aan het eind aangeeft dat het naar het volgende deeltje gaat, moet die code alleen stoppen bij het laatste deeltje.
We hebben hier gebruik gemaakt van het `if` statement.
Deze manier van werken is misschien niet heel efficient - het kost rekentijd - maar kan enorm helpen bij debuggen van code: je ziet snel welke tekst wel en welke tekst niet geprint wordt.



```{exercise}
:label: intro-ex-2
Maak 10 deeltjes, elk startend op (een andere positie op) de x-as. Elk deeltje heeft een eigen random snelheid tussen de -5 en +5.
Update de positie en plot de positie van elk deeltje. 

En als je wilt gaan voor een uitdaging: 
- plot in dezelfde figuur de begin en eindpositie van elk deeltje
- maak gebruik van een kleine verticale snelheid om het verschil te kunnen zien.
- hoe zou je elk deeltje een andere kleur kunnen geven in de plot zodat ze ook echt traceerbaar zijn?
```

```python
import random
x_waardes = np.arange(10)
class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m
        self.v = np.array(v, dtype=float)
        self.r = np.array(r, dtype=float)
        self.R = np.array(R, dtype=float)

    def update_position(self, dt):
        self.r += self.v * dt

particles = [ParticleClass(1, [random.uniform(-5,5), random.uniform(-1,1)], [x_waardes[i], 0], 1) for i in range(10)]
colors = plt.cm.tab10(range(10))
for i, p in enumerate(particles):
    start = p.r.copy()      
    p.update_position(1.0) 
    
    plt.plot([start[0], p.r[0]], [start[1], p.r[1]], '-', color=colors[i])
    plt.scatter(*start, color=colors[i], marker='o') 
    plt.scatter(*p.r, color=colors[i], marker='x')  
plt.show()
```

Onze simulaties worden snel “zwaar”. Er moeten veel berekeningen gedaan worden en dat kost nu eenmaal tijd. Nu is het wel zo dat sommige berekeningen meer tijd kosten dan anderen. Om te kijken of een stukje code geoptimaliseerd kan worden, moeten we een timer maken. Dat kan met onderstaande code:

```python
import time
import numpy as np

N = int(1e5)          # aantal random getallen
herhalingen = 5       # aantal keer dat we meten, voor een betere schatting

# Methode 1: element voor element toevoegen (traag: elke append kopieert het hele array)
tijden1 = np.zeros(herhalingen)
for k in range(herhalingen):
    start = time.perf_counter()
    array = np.array([])
    for i in range(N):
        array = np.append(array, np.random.rand(1))
    eind = time.perf_counter()
    tijden1[k] = eind - start
len1 = np.mean(tijden1)
print("Methode 1 kostte", len1, "seconde!")

# Methode 2: in één keer N random getallen maken
tijden2 = np.zeros(herhalingen)
for k in range(herhalingen):
    start = time.perf_counter()
    array = np.random.rand(N)
    eind = time.perf_counter()
    tijden2[k] = eind - start
len2 = np.mean(tijden2)
print("Methode 2 kostte", len2, "seconde!")

keer_korter = float(len1/len2)
print("Methode 2 is", round(keer_korter), "keer sneller")
```

Idealiter draaien we de bovenstaande code een aantal keer zodat we een goede schatting hebben van de rekentijd.


```{exercise} Timen
:label: ex-intro-3

Bovenstaande code kan natuurlijk veel sneller, in de functie random.rand() kun je het aantal elementen (N) mee geven. D
Voeg in bovenstaande code-cell de alternatieve methode voor het maken van de array toe en bereken de runtime.
Hoeveel keer korter is de runtime?
```

```{tip}
Met de MyST extentie in jupyterlab kun je gebruik maken van [eval](https://mystmd.org/guide/notebooks-with-markdown#myst-inline-expressions). Hieronder kun je de naam van je variabele invoeren die je hierboven heb gespecificeerd om uit te drukken hoeveel korter de runtime is.
```

```{solution} ex-intro-3
De runtime is {eval}`keer_korter` korter.
```


## Aan de slag met deeltjes
Zoom je heel ver in, dan zie je deeltjes rond vliegen.
Elk met een eigen massa, een eigen snelheid en richting.
De deeltjes botsen onderling, wisselen energie uit.
We zouden op basis van een botsingsmodel van deeltjes iets moeten kunnen leren over thermodynamica.
In de thermodynamica gaat het dan om heel veel deeltjes. Maar laten we beginnen met twee botsende 'deeltjes'.

```{exercise}
:label: ex-deeltjes-1
Bekijk onderstaande video.
```

```{iframe} https://www.youtube.com/embed/HEfHFsfGXjs?si=g37KWYOi0LqmuAkD
```

Je raadt het misschien al... deze simulatie gaan we na bouwen!
Daarbij maken we gebruik van Python classes uit het vorige hoofdstuk en de basics van simulaties geleerd in Q1.
Zorg er dus voor dat je weet hoe dit werkt!
We maken ook gebruik van plotten, en daarbovenop een animatie.
Hoe de animatie precies werkt en hoe je die zelf maakt hoef je niet te weten.
Je zou wel in staat moeten zijn de code te lezen.

In de onderstaande cell maken we de ParticleClass aan, en geven we enkele parameters van onze simulatie op.

```python
# Importeren van libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Maken van de class
class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m                         
        self.v = np.array(v, dtype=float)  
        self.r = np.array(r, dtype=float)  
        self.R = R

    def update_position(self):
        self.r += self.v * dt

# Simulation parameters
dt = 0.1                                                            # tijd stap
num_steps = 500                                                     # aantal te nemen stappen
particle = ParticleClass(m=1.0, v=[5.0, 0], r=[0.0, 0.0], R=1.0)    # het maken van ons deeltje
```

```{note} dt
In bovenstaande `update_position` functie hebben we de tijdstap niet als inputvariabele meegegeven, waarbij we dat in de eerdere module wel hebben gedaan.
In deze module gaan we er vanuit dat de tijdstap altijd hetzelfde is voor alle deeltjes in de simulatie.
Daarom definiëren we de tijdstap als een globale variabele in de simulatie.
```

We hebben nu een deeltje met massa, een snelheid, een begin positie en een straal.
We hebben ook al de stapgrootte bepaald!

We willen de beweging van dat deeltje straks bestuderen en moeten dus een plot maken:

```python
# Creeer een figuur en de assen
fig, ax = plt.subplots()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Animatie")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Toon het deeltje als een rode stip
dot, = ax.plot([], [], 'ro', markersize=10); # semicolon to suppress output
```

Bij het aanmaken van ons deeltje hebben we het deeltje een beginpositie en snelheid mee gegeven.
Als we dan per tijdstap de positie bepalen en deze laten plotten en die plots achter elkaar plakken, dan krijgen we een animatie van het deeltje.
Met FuncAnimation wordt die animatie voor ons gedaan.

```python
# Initialisieren van de functie voor de animatie
def init():
    dot.set_data([], [])
    return dot,

# Updaten van de functie voor elk frame
def update(frame):
    particle.update_position()
    dot.set_data([particle.r[0]], [particle.r[1]])
    return dot,

# Creeer de animatie
ani = FuncAnimation(fig, update, frames=range(200), init_func=init, blit=True, interval=50)

# Omdat we werken met Jup. Notebooks (en niet een .py file)
from IPython.display import HTML
HTML(ani.to_jshtml())

```

`Frames` in FuncAnimation verwijst naar het totaal aantal frames dat gebruikt wordt.
`Interval` naar de snelheid van de animatie, nl. 50 milliseconden ofwel 20 frames per seconde.

Merk op dat als we de laatste cel opnieuw runnen, het deeltje zich niet in de oorsprong bevindt.
Dat is een 'eigenaardigheid' van Jupyter Notebooks.
Het is nu beter om alle code in één cel te plaatsen (hieronder gedaan voor je), zodat we ervoor zorgen dat het deeltje altijd in de oorsprong begint.


```python
# Simulation parameters
dt = 0.1         # time step
num_steps = 500  # number of time steps
particle = ParticleClass(m=1.0, v=[5.0, 0], r=[0.0, 0.0],R=1.0)  

# Create the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Create the particle as a red dot
dot, = ax.plot([], [], 'ro', markersize=10)

# Initialization function for animation
def init():
    dot.set_data([], [])
    return dot,

# Update function for each frame
def update(frame):
    particle.update_position()
    dot.set_data([particle.r[0]], [particle.r[1]])
    return dot,

# Create animation
ani = FuncAnimation(fig, update, frames=range(200), init_func=init, blit=True, interval=50)

# For Jupyter notebook:
from IPython.display import HTML
HTML(ani.to_jshtml())

```

Een van de dingen die je kunt opmerken, is dat het deeltje niet in zijn doos blijft.

```{exercise} Doorlopende doos
:label: ex-deeltjesmodel-2
Pas de onderstaande code aan zodat de doos een aaneengesloten doos is; als hij links eruit vliegt, komt hij er rechts in (en vice versa).

Controleer je eigen antwoord door de functie hierboven even te vervangen door deze functie, bekijk het resultaat.
```

```python tags=["NB1_doorlopendedoos"]
# doorlopende doos
def update(frame):
    particle.update_position()
    dot.set_data([particle.r[0]], [particle.r[1]])

    if particle.r[0] > 10:
        particle.r[0] = -10
    elif particle.r[0] < -10:
        particle.r[0] = 10

    if particle.r[1] > 10:
        particle.r[1] = -10
    elif particle.r[1] < -10:
        particle.r[1] = 10
### YOUR CODE HERE ###

    return dot,         
```

```{exercise} Harde wanden
:label: ex-deeltjesmodel-3
Een tweede optie is dat we een doos hebben met harde wanden. Op het moment dat het deeltje de wand raakt, wordt deze gereflecteerd. Schrijf hieronder de code zodat het deeltje in zijn doos blijft, waarbij de doos harde wanden heeft.

Test wederom je code door de vernieuwde functie hierboven te plakken.
```

```python tags=["NB1_hardewand"]
# harde wanden
def update(frame):
    particle.update_position()
    dot.set_data([particle.r[0]], [particle.r[1]])
    # Alleen omkeren als het deeltje de wand raakt én naar buiten beweegt.
    # Zo blijft het deeltje niet 'hangen' in de wand als het in één stap net voorbij de wand schiet.
    if (particle.r[0] >= 10 and particle.v[0] > 0) or (particle.r[0] <= -10 and particle.v[0] < 0):
        particle.v[0] *= -1

    if (particle.r[1] >= 10 and particle.v[1] > 0) or (particle.r[1] <= -10 and particle.v[1] < 0):
        particle.v[1] *= -1
    return dot,
```

```{exercise}
:label: ex-deeltjesmodel-4

Verander nu de code zodanig dat de beginsnelheid gegeven wordt door $\vec{v}=5\hat{x} + 3\hat{y}$ en het deeltje botst tegen alle muren.
```

```python
# Simulation parameters
dt = 0.1         # time step
num_steps = 500  # number of time steps

particle = ParticleClass(m=1.0, v=[5.0, 3.0], r=[0.0, 0.0], R=1.0)   # v = 5 x + 3 y

# Create the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Create the particle as a red dot
dot, = ax.plot([], [], 'ro', markersize=10);

# Initialization function for animation
def init():
    dot.set_data([], [])
    return dot,

# Update function for each frame
def update(frame):
    particle.update_position()
    dot.set_data([particle.r[0]], [particle.r[1]])
    # harde wanden: omkeren als het deeltje de wand raakt en naar buiten beweegt
    for j in range(2):   # j = 0: x-richting, j = 1: y-richting
        if (particle.r[j] >= 10 and particle.v[j] > 0) or (particle.r[j] <= -10 and particle.v[j] < 0):
            particle.v[j] *= -1
    return dot,

# Create animation
ani = FuncAnimation(fig, update, frames=range(200), init_func=init, blit=True, interval=50)

# For Jupyter notebook:
from IPython.display import HTML
HTML(ani.to_jshtml())
```

Laten we teruggaan naar ons deeltje. Er is een functie om de positie bij te werken, hoewel de snelheid hetzelfde lijkt te blijven... kunnen we de snelheid veranderen door (bijvoorbeeld) de versnelling door de zwaartekracht?

```{exercise}
:label: ex-deeltjesmodel-5

Hieronder is de ParticleClass aangepast zodat er gebruik gemaakt kan worden van een versnelling.
Maak de simulatie van het deeltje zodat het zich beweegt in een zwaartekrachtsveld met $a = -9.81\hat{y}$.

*Tip:* De tweede solution kun je van hierboven kopieren natuurlijk!
```

```python
# Maken van de class met versnelling
class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m                         # mass of the particle
        self.v = np.array(v, dtype=float)  # velocity vector
        self.r = np.array(r, dtype=float)  # position vector
        self.R = R                         # radius of the particle

    def update_position(self):
        self.r += self.v * dt
    
    def update_velocity(self, a):
        self.v += a*dt
```


```python

# Simulation parameters
dt = 0.1         # time step
num_steps = 500  # number of time steps
particle = ParticleClass(m=1.0, v=[5.0, 0], r=[0.0, 0.0],R=1.0)  
a = np.array([0.0, -9.81])   # zwaartekrachtsversnelling (m/s^2)
# Create the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Create the particle as a red dot
dot, = ax.plot([], [], 'ro', markersize=10);

# Initialization function for animation
def init():
    dot.set_data([], [])
    return dot,

# Update function for each frame
def update(frame):
    particle.update_velocity(a)
    particle.update_position()
    
    dot.set_data([particle.r[0]], [particle.r[1]])
    # harde wanden
    for j in range(2):
        if (particle.r[j] >= 10 and particle.v[j] > 0) or (particle.r[j] <= -10 and particle.v[j] < 0):
            particle.v[j] *= -1
    return dot,

# Create animation
ani = FuncAnimation(fig, update, frames=range(200), init_func=init, blit=True, interval=50)

# For Jupyter notebook:
from IPython.display import HTML
HTML(ani.to_jshtml())

```

```{exercise}
:label: ex-deeltjesmodel-6

Wat valt je op aan het behoud van energie?!
Herinner wat we in Q1 hebben gedaan bij de basics van een numerieke simulatie...
```

Om dit te controleren rekenen we dezelfde simulatie door zonder animatie en houden we de totale energie per massa bij: $E/m = \tfrac{1}{2}|\vec v|^2 + g\,y$.

```python
def energie_simulatie(dt, t_eind=20):
    deeltje = ParticleClass(m=1.0, v=[5.0, 0], r=[0.0, 0.0], R=1.0)
    a = np.array([0.0, -9.81])
    n = int(t_eind/dt)
    t = np.arange(n)*dt
    E = np.zeros(n)
    for i in range(n):
        deeltje.v += a*dt          # zelfde volgorde als in de animatie
        deeltje.r += deeltje.v*dt
        for j in range(2):
            if (deeltje.r[j] >= 10 and deeltje.v[j] > 0) or (deeltje.r[j] <= -10 and deeltje.v[j] < 0):
                deeltje.v[j] *= -1
        E[i] = 0.5*np.dot(deeltje.v, deeltje.v) - a[1]*deeltje.r[1]   # energie per kg
    return t, E

plt.figure()
for dt_test in [0.1, 0.01, 0.001]:
    t_E, E = energie_simulatie(dt_test)
    plt.plot(t_E, E, label=f'$\\Delta t$ = {dt_test} s')
plt.axhline(12.5, color='k', ls='--', label='exacte beginenergie')
plt.xlabel('$t$ (s)')
plt.ylabel('$E/m$ (J/kg)')
plt.legend()
plt.show()
```

```{solution} ex-deeltjesmodel-6
De energie is **niet** behouden, terwijl dat natuurkundig wel zou moeten (er werkt alleen de zwaartekracht en de botsingen met de wanden zijn elastisch). Het deeltje stuitert daardoor niet steeds even hoog: bij elke stuiter verandert de energie een stukje (in deze code wordt de stuiterhoogte steeds lager, met een grotere $\Delta t$ gaat dat sneller).

Dit komt door de numerieke methode (Euler), net als bij de fietsers in Q1: binnen een tijdstap houden we de snelheid constant, terwijl die in werkelijkheid continu verandert. Daarnaast schiet het deeltje in één stap een stukje voorbij de wand ($|v|\Delta t$) voordat de snelheid wordt omgedraaid, waardoor het op de verkeerde hoogte wordt teruggekaatst. Beide fouten zijn evenredig met $\Delta t$: hoe kleiner de tijdstap, hoe beter de energie behouden blijft (zie de figuur).
```

Een optie om de simulatie te verbeteren is de tijdstap $\Delta t$ kleiner te maken, maar dan hebben we ook meer geduld meer nodig - het aantal berekeningen schaalt met $\frac{1}{\Delta t}$. Een tweede optie is een meer directe oplossing: We weten dat $a = \mathrm{const.}$ en daarom weten we ook de bewegingsvergelijking van het deeltje!

We willen ook graag weten waar het deeltje is geweest, dat is in onderstaande code toegevoegd.

```python
# Maken van de class met versnelling
class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m                  
        self.v = np.array(v, dtype=float)  
        self.r = np.array(r, dtype=float)  
        self.R = R  

    def update_position(self):
        self.r += self.v * dt + 1/2 * a * dt**2  
    
    def update_velocity(self, a):
        """Update the particle's velocity."""
        self.v += a*dt

# Simulation parameters
dt = 0.1         
num_steps = 500  
particle = ParticleClass(m=1.0, v=[0, 0], r=[0.0, 0.0],R=1.0)  
a = np.array([0.0, -5.0])  

track_x = []
track_y = []

# creeeren van de plot en de assen
fig, ax = plt.subplots()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")
track_line, = ax.plot([], [], 'r--', linewidth=1)  

# creeeren van ons rode deeltje
dot, = ax.plot([], [], 'ro', markersize=10);

# initializeren van onze functie voor de animatie
def init():
    dot.set_data([], [])
    return dot,

# Update function for each frame
def update(frame):
    particle.update_position()
    particle.update_velocity(a)

    track_x.append(particle.r[0])
    track_y.append(particle.r[1])
    track_line.set_data(track_x, track_y)
    
    dot.set_data([particle.r[0]], [particle.r[1]])
    if particle.r[0]**2>100: # Check if particle is outside the bounds, np.abs could be used but is slower
        particle.v[0] = -particle.v[0]
    if particle.r[1]**2>100: # Check if particle is outside the bounds, np.abs could be used but is slower
        particle.v[1] = -particle.v[1]
    return dot, track_line

# Create animation
ani = FuncAnimation(fig, update, frames=range(200), init_func=init, blit=True, interval=50)

# For Jupyter notebook:
from IPython.display import HTML
HTML(ani.to_jshtml())
```

We hebben steeds slechts gewerkt met een enkel deeltje. Maar om de simulatie uit het filmpje te maken, hebben we twee deeltjes nodig.

```{exercise}
:label: ex-deeltjesmodel-7

Maak de code hieronder voor twee deeltjes die bewegen in een afgesloten doos. Beide deeltjes ondervinden een zwaartekracht. Beide deeltjes starten op hetzelfde punt, zonder beginsnelheid in de verticale richting. Echter, een van de twee deeltjes heeft een initiële horizontale snelheid.
```

```python
# Simulation parameters
dt = 0.1         # time step
num_steps = 500  # number of time steps
particleA = ParticleClass(m=1.0, v=[0, 0], r=[0.0, 0.0],R=1.0)  
particleB = ParticleClass(m=1.0, v=[4.0, 0], r=[0.0, 0.0],R=1.0)   # zelfde startpunt, wel een horizontale beginsnelheid
a = np.array([0.0, -9.81])                                          # zwaartekracht werkt op beide deeltjes

track_x = []
track_y = []
trackB_x = []
trackB_y = []

# Create the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")
track_line, = ax.plot([], [], 'r--', linewidth=1)  
trackB_line, = ax.plot([], [], 'b--', linewidth=1)  

# Create the particle as a red dot
dotA, = ax.plot([], [], 'ro', markersize=10)
dotB, = ax.plot([], [], 'bo', markersize=10)

# Initialization function for animation
def init():
    dotA.set_data([], [])
    dotB.set_data([], [])
    return dotA, dotB

# Update function for each frame
def update(frame):
    particleA.update_position()
    particleA.update_velocity(a)

    particleB.update_position()
    particleB.update_velocity(a)

    track_x.append(particleA.r[0])
    track_y.append(particleA.r[1])
    track_line.set_data(track_x, track_y)

    trackB_x.append(particleB.r[0])
    trackB_y.append(particleB.r[1])
    trackB_line.set_data(trackB_x, trackB_y)
    
    dotA.set_data([particleA.r[0]], [particleA.r[1]])
    dotB.set_data([particleB.r[0]], [particleB.r[1]])

    # harde wanden voor beide deeltjes: alleen omkeren als het deeltje naar buiten beweegt
    for p in [particleA, particleB]:
        for j in range(2):
            if (p.r[j] >= 10 and p.v[j] > 0) or (p.r[j] <= -10 and p.v[j] < 0):
                p.v[j] = -p.v[j]

    return dotA, dotB, track_line, trackB_line

# Create animation
ani = FuncAnimation(fig, update, frames=range(200), init_func=init, blit=True, interval=50)

# For Jupyter notebook:
from IPython.display import HTML
HTML(ani.to_jshtml())
```

Nu je een beetje weet hoe een numeriek model werkt, hoe je boundaries kunt opgeven, deeltjes kunt laten versnellen etc, gaan we de deeltjes echt laten interacteren. Hierboven bewogen ze immers nog door elkaar heen...

## Botsingen in 1D
```{exercise}
:label: ex-deeltjes_10

1. Gegeven twee particles elk met straal $R$ en posities $r_1$ en $r_2$. Wat is de (wiskundige) voorwaarde voor een botsing?

2. Kijk naar de onderstaande code en schrijf de voorwaarde voor een botsing op (m.a.w. maak de functie `collide_detection(self, other)` werkend).

3. Vul onderstaande code aan waarbij de twee deeltjes van richting veranderen wanneer ze botsen.

4. Automatiseer en maak de code eleganter waar mogelijk. Kijk daarvoor ook terug naar het aanmaken van een particlearray.
```

```{solution} ex-deeltjes_10
1. Twee bolletjes raken elkaar als de afstand tussen hun middelpunten kleiner of gelijk is aan de som van de stralen:
$$|\vec r_1 - \vec r_2| \le R_1 + R_2 \quad (=2R \text{ bij gelijke stralen}).$$
Om de wortel te vermijden (scheelt rekentijd) kun je ook beide kanten kwadrateren: $(x_1-x_2)^2 + (y_1-y_2)^2 \le (R_1+R_2)^2$.

2. Zie `collide_detection` hieronder: `dx`, `dy` zijn de verschillen in positie en `rr` de som van de stralen.

3. Bij een botsing keren beide deeltjes hun snelheid om. Dat doen we alleen als de deeltjes naar elkaar toe bewegen, $(\vec v_1-\vec v_2)\cdot(\vec r_1-\vec r_2) < 0$. Anders kan het gebeuren dat de deeltjes na de botsing in de volgende tijdstap nog steeds overlappen, waardoor ze opnieuw 'botsen' en aan elkaar blijven plakken.

4. In de elegantere versie (cel daaronder) staan de deeltjes in een lijst, wordt de wandbotsing een methode van de class, en lopen we met een dubbele `for`-loop over alle paren deeltjes. Zo werkt dezelfde code voor 2 of voor 100 deeltjes.
```

```python tags=["NB1_botsingsvoorwaarde"]
# Define a class for a particle

# Maken van de class met botsing
class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m                  
        self.v = np.array(v, dtype=float)  
        self.r = np.array(r, dtype=float)  
        self.R = R 

    def update_position(self):
        self.r += self.v * dt  

    def collide_detection(self, other):
        dx = self.r[0] - other.r[0]    # verschil in x-positie
        dy = self.r[1] - other.r[1]    # verschil in y-positie
        rr = self.R + other.R          # som van de stralen
        return  dx**2+dy**2 < rr**2 
```

```python
# Simulation parameters
dt = 0.1         # time step
num_steps = 200  # number of time steps

particleA = ParticleClass(m=1.0, v=[2.5, 0], r=[-2.0, 0.0],R=0.45)  
particleB = ParticleClass(m=1.0, v=[-1, 0], r=[0.0, 0.0],R=0.45)  

track_x = []
track_y = []


# Creer de plot

fig, ax = plt.subplots()

ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")
track_line, = ax.plot([], [], 'r--', linewidth=1)  

# Toon het deeltje als een rode stip
dot, = ax.plot([], [], 'ro', markersize=10); # semicolon to suppress output
dotA, = ax.plot([], [], 'ro', markersize=10)
dotB, = ax.plot([], [], 'bo', markersize=10)

# Initaliseren voor de animatie
def init():
    dot.set_data([], [])
    return dot,

# Updaten van de functie per frame
def update(frame):
    particleA.update_position()
    particleB.update_position()

    track_x.append(particleA.r[0])
    track_y.append(particleA.r[1])
    track_line.set_data(track_x, track_y)
    
    dotA.set_data([particleA.r[0]], [particleA.r[1]])
    dotB.set_data([particleB.r[0]], [particleB.r[1]])

    # botsing tussen de deeltjes onderling
    if particleA.collide_detection(particleB):
        # alleen omkeren als de deeltjes naar elkaar toe bewegen (voorkomt 'plakken')
        if np.dot(particleA.v - particleB.v, particleA.r - particleB.r) < 0:
            particleA.v = -particleA.v
            particleB.v = -particleB.v

    # botsing met de wand
    if particleA.r[0]**2>100: 
        particleA.v[0] = -particleA.v[0]
    if particleA.r[1]**2>100: 
        particleA.v[1] = -particleA.v[1]

    dot.set_data([particleB.r[0]], [particleB.r[1]])
    if particleB.r[0]**2>100: 
        particleB.v[0] = -particleB.v[0]
    if particleB.r[1]**2>100: 
        particleB.v[1] = -particleB.v[1]

    return dot, track_line

# Creeer animatie
ani = FuncAnimation(fig, update, frames=range(num_steps), init_func=init, blit=True, interval=50)

# Voor Jupyter notebook:
from IPython.display import HTML
HTML(ani.to_jshtml())
```

Elegantere versie (onderdeel 4): de deeltjes staan in een lijst en de wandbotsing is een methode van de class. Dezelfde code werkt nu voor een willekeurig aantal deeltjes.

```python
class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m
        self.v = np.array(v, dtype=float)
        self.r = np.array(r, dtype=float)
        self.R = R

    def update_position(self):
        self.r += self.v * dt

    def collide_detection(self, other):
        return np.sum((self.r - other.r)**2) < (self.R + other.R)**2

    def approaching(self, other):
        """True als de deeltjes naar elkaar toe bewegen."""
        return np.dot(self.v - other.v, self.r - other.r) < 0

    def wall_collision(self, L=10):
        """Harde wanden op -L en +L in zowel x als y."""
        naar_buiten = ((self.r >= L - self.R) & (self.v > 0)) | ((self.r <= -L + self.R) & (self.v < 0))
        self.v[naar_buiten] *= -1

# Simulatie parameters
dt = 0.1
num_steps = 200
kleuren = ['r', 'b']

particles = [ParticleClass(m=1.0, v=[2.5, 0], r=[-2.0, 0.0], R=0.45),
             ParticleClass(m=1.0, v=[-1, 0],  r=[0.0, 0.0],  R=0.45)]

fig, ax = plt.subplots()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")
dots = [ax.plot([], [], kleuren[k % len(kleuren)] + 'o', markersize=10)[0] for k in range(len(particles))]

def init():
    for d in dots:
        d.set_data([], [])
    return dots

def update(frame):
    for p in particles:
        p.update_position()
        p.wall_collision()

    # botsingen tussen alle paren deeltjes (elk paar één keer)
    for i in range(len(particles)):
        for j in range(i + 1, len(particles)):
            pi, pj = particles[i], particles[j]
            if pi.collide_detection(pj) and pi.approaching(pj):
                pi.v = -pi.v
                pj.v = -pj.v

    for p, d in zip(particles, dots):
        d.set_data([p.r[0]], [p.r[1]])
    return dots

ani = FuncAnimation(fig, update, frames=range(num_steps), init_func=init, blit=True, interval=50)

from IPython.display import HTML
HTML(ani.to_jshtml())
```

Terug naar ons vraagstuk... 
we willen een simulatie waarbij een deeltje met een massa $m_1$ en snelheid $v_1$ op een andere stilstaand deeltje met massa $m_2$ botst.
Deeltje twee beweegt naar een muur, botst tegen de muur en beweegt richting deeltje 1 en botst tegen dit deeltje.
Hoe vaak vindt deze botsing plaats als functie van de massa verhouding $\frac{m_1}{m_2}$?

Daarvoor moeten we even terug naar het botsingsmodel zoals geleerd in [*Klassieke Mechanica*.](https://freekpols.github.io/Mechanica/content/classic/collisions)
Bij elastische botsingen is zowel het impulsmoment als de kinetische energie behouden. 

$$
    \sum_i \vec{p}_i^{before} = \sum_i \vec{p}_i^{after}
$$

$$ 
    E_{kin, before} = E_{kin, after} 
$$

Voor een botsing met twee deeltjes levert dit een analytische oplossing:

$$
\begin{align}
\mathbf{v}'_1 &= \mathbf{v}_1-\frac{2 m_2}{m_1+m_2} \ \frac{\langle \mathbf{v}_1-\mathbf{v}_2,\,\mathbf{x}_1-\mathbf{x}_2\rangle}{\|\mathbf{x}_1-\mathbf{x}_2\|^2} \ (\mathbf{x}_1-\mathbf{x}_2),
\\
\mathbf{v}'_2 &= \mathbf{v}_2-\frac{2 m_1}{m_1+m_2} \ \frac{\langle \mathbf{v}_2-\mathbf{v}_1,\,\mathbf{x}_2-\mathbf{x}_1\rangle}{\|\mathbf{x}_2-\mathbf{x}_1\|^2} \ (\mathbf{x}_2-\mathbf{x}_1)
\end{align}
$$

we vragen je natuurlijk niet om deze vergelijking zelf te schrijven in Python, maar die vergelijking operationaliseren we hieronder wel.

```{exercise} Lezen van code
:label: ex-deeltjes-10

Lees onderstaande code, en probeer wijs te worden wat elk deel van de code doet. 
In het volgende hoofdstuk wordt de code net iets anders gestructureerd en uitgebreider uitgelegd.
```

```{exercise} Toevoegen van een counter
:label: ex-deeltjes-11

Als je de code runt, zie je dat het botsingsmodel goed lijkt te werken.
Er is zelfs een counter toegevoegd die het aantal botsingen dat plaats vindt zou moeten weergeven.
Maar in de code wordt de counter nog niet geupdate.

Voeg op de juiste plaats(en) een update van de counter toe (per botsing moet de counter 1 omhoog gezet worden).
```

```{exercise} Controleer
:label: ex-deeltjes-12

Controleer de uitspraak in de video dat $\pi$ in eens opduikt bij de massaverhouding. 

Wat is de beperkende factor in ons script om het aantal botsingen als functie van de massaverhouding te onderzokeen?
```

```{solution} ex-deeltjes-12
Met de code onderaan (`tel_botsingen`) vinden we:

| $m_2/m_1$ | 1 | 100 | 10 000 |
|---|---|---|---|
| aantal botsingen | 3 | 31 | 314 |

Voor $m_2/m_1 = 100^n$ geeft het aantal botsingen precies de eerste $n+1$ cijfers van $\pi$, zoals in de video. (Dit volgt uit behoud van energie: in de coördinaten $(\sqrt{m_1}v_1, \sqrt{m_2}v_2)$ ligt de toestand op een cirkel en elke botsing draait de toestand over een vaste hoek $2\arctan\sqrt{m_1/m_2}$.)

De beperkende factor is de **tijdstap** $\Delta t$ (en daarmee de rekentijd). Bij een grote massaverhouding volgen de botsingen tussen het lichte deeltje, het zware deeltje en de muur elkaar heel snel op. Als het lichte deeltje in één tijdstap verder beweegt dan de ruimte tussen muur en zwaar deeltje, worden botsingen gemist of te laat gedetecteerd en klopt de telling niet meer. Daarom moet $\Delta t$ kleiner bij een grotere massaverhouding (voor $10\,000$ was $\Delta t = 0{,}001$ nodig), en schaalt het aantal rekenstappen met $1/\Delta t$. Bij nog grotere verhoudingen speelt ook de afrondingsfout van floating point getallen mee.
```


```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class ParticleClass:
    def __init__(self, m, v, r, R):
        self.m = m                  # mass of the particle
        self.v = np.array(v, dtype=float)  # velocity vector
        self.r = np.array(r, dtype=float)  # position vector
        self.R = R                         # radius of the particle

    def update_position(self):
        self.r += self.v * dt 

    def collide_detection(self, other):
        return np.linalg.norm(self.r - other.r) <= (self.R + other.R)

# Simulation parameters
dt = 0.1         # time step
num_steps = 300  # number of time steps
m1 = 1.0
m2 = 50.0        # opdracht hieronder: m1 = 1 en m2 = 50

particleA = ParticleClass(m=m1, v=[0, 0], r=[-4.0, 0.0],R=0.45)  
particleB = ParticleClass(m=m2, v=[-1, 0], r=[-2.0, 0.0],R=0.45)  


# Create the figure and axis

fig, ax = plt.subplots()

ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_aspect('equal')
ax.set_title("Particle Animation")
ax.set_xlabel("x")
ax.set_ylabel("y")
dot, = ax.plot([], [], 'ro', markersize=10); # semicolon to suppress output

counter = 0

# Create the particle as a red dot
dotA, = ax.plot([], [], 'ro', markersize=10)
dotB, = ax.plot([], [], 'bo', markersize=10)

counter_text = ax.text(-9.5, 9, "")

# Initialization function for animation
def init():
    dot.set_data([], [])
    return dot,

# Update function for each frame
def update(frame):
    global counter
    particleA.update_position()
    particleB.update_position()

    dotA.set_data([particleA.r[0]], [particleA.r[1]])
    dotB.set_data([particleB.r[0]], [particleB.r[1]])

    #collision detection and response
    # alleen botsen als de deeltjes naar elkaar toe bewegen, anders telt één botsing
    # meerdere keren mee zolang de deeltjes nog overlappen
    if particleA.collide_detection(particleB) and np.dot(particleA.v - particleB.v, particleA.r - particleB.r) < 0:
        vA, vB, mA, mB, rA, rB = particleA.v, particleB.v, particleA.m, particleB.m, particleA.r, particleB.r
        vA_new = vA - 2 * mB / (mA + mB) * np.dot(vA - vB, rA - rB) / (1e-12+np.linalg.norm(rA - rB))**2 * (rA - rB)
        vB_new = vB - 2 * mA / (mA + mB) * np.dot(vB - vA, rB - rA) / (1e-12+np.linalg.norm(rB - rA))**2 * (rB - rA)
        particleA.v = vA_new
        particleB.v = vB_new
        counter += 1


    # wall collision detection and response
    # Zoals in de video is er alleen een muur aan de linkerkant (x = -10).
    # Alleen omkeren als het deeltje de muur raakt én naar de muur toe beweegt.
    if particleA.r[0] < -10 and particleA.v[0] < 0:
        particleA.v[0] = -particleA.v[0]
        counter += 1
    if particleA.r[1]**2>100 and particleA.r[1]*particleA.v[1] > 0: 
        particleA.v[1] = -particleA.v[1]
        counter += 1

    if particleB.r[0] < -10 and particleB.v[0] < 0: 
        particleB.v[0] = -particleB.v[0]
        counter += 1
    if particleB.r[1]**2>100 and particleB.r[1]*particleB.v[1] > 0: 
        particleB.v[1] = -particleB.v[1]
        counter += 1

    dot.set_data([particleB.r[0]], [particleB.r[1]])
    counter_text.set_text(f"Collisions: {counter}")

    return dotA, dotB, counter_text

# Create animation
ani = FuncAnimation(fig, update, frames=range(num_steps), init_func=init, blit=True, interval=50)

# For Jupyter notebook:
from IPython.display import HTML
HTML(ani.to_jshtml())
```

```{exercise}
Pas de bovenstaande code aan zodat $m_1 = 1$ en $m_2 = 50$. Run de code en geef hieronder als integer aan hoeveel botsingen totaal hebben plaats gevonden. Let op! Het kan zijn dat je voor het precies antwoord nog wat dingen moet aanpassen...
```


De animatie hierboven (met $m_2 = 50$) eindigt op 22 botsingen. Om het antwoord niet af te laten hangen van het aantal frames, rekenen we hetzelfde nog eens door zonder animatie, net zo lang tot er geen botsingen meer kunnen komen: beide deeltjes bewegen dan naar rechts en het lichte deeltje A is langzamer dan B.

```python
def tel_botsingen(m1, m2, dt=0.01):
    """Telt het aantal botsingen (deeltje-deeltje + deeltje-muur) voor de opstelling hierboven."""
    A = ParticleClass(m=m1, v=[0, 0],  r=[-4.0, 0.0], R=0.45)
    B = ParticleClass(m=m2, v=[-1, 0], r=[-2.0, 0.0], R=0.45)
    n = 0
    while True:
        A.r += A.v*dt
        B.r += B.v*dt
        if A.collide_detection(B) and np.dot(A.v - B.v, A.r - B.r) < 0:
            vA = A.v - 2*B.m/(A.m + B.m)*np.dot(A.v - B.v, A.r - B.r)/np.linalg.norm(A.r - B.r)**2*(A.r - B.r)
            vB = B.v - 2*A.m/(A.m + B.m)*np.dot(B.v - A.v, B.r - A.r)/np.linalg.norm(B.r - A.r)**2*(B.r - A.r)
            A.v, B.v = vA, vB
            n += 1
        if A.r[0] < -10 and A.v[0] < 0:     # muur links
            A.v[0] = -A.v[0]
            n += 1
        # klaar: beide naar rechts en A haalt B nooit meer in
        if A.v[0] >= 0 and B.v[0] >= 0 and A.v[0] <= B.v[0]:
            return n

for m2_test, dt_test in [(1, 0.01), (100, 0.01), (10000, 0.001)]:
    print(f"m2/m1 = {m2_test:>6}: {tel_botsingen(1.0, m2_test, dt_test)} botsingen")
```

```python tags=["NB1_eindantwoord"]
num_botsingen = tel_botsingen(1.0, 50.0)
print(num_botsingen)     # 22
```

```{exercise}
Klaar? Bekijk de twee versies, kies de beste en doe een `Pull Request` naar `Main`. Je partner reviewt de code, doet eventueel suggesties. Je partner is ook degene die uiteindelijke  de Pull Request goed keurt.
```
