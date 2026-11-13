import networkx as nx
from Bio.PDB import *
import warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mayavi.mlab as mlab
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
from tkinter import simpledialog
from tkinter import filedialog, messagebox
import os
import requests


from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def BuscarArchivo():
    global filename
    filename = filedialog.askopenfilename(initialdir="/", title="Select file", filetypes=(("Archivos pdb", "*.pdb"), ("Todos los archivos", "*.*")))
    if filename:
        label_archivo.config(text=f"Selected file: {filename}")
    
    filtered_lines = []

    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith("HETATM") and not line.startswith("CONECT") and not line.startswith("SOL"):
                filtered_lines.append(line)

    with open(filename, 'w') as f:
        f.writelines(filtered_lines)
        

def descargar_pdb_directo(pdb_id, destino=None):
    """
    Descarga un archivo PDB directamente desde el servidor RCSB PDB.

    Args:
        pdb_id (str): Identificador del archivo PDB.
        destino (str, opcional): Directorio donde guardar el archivo.

    Returns:
        str: Ruta del archivo descargado, o None si ocurre un error.
    """
    base_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    destino = destino or os.getcwd()
    archivo_destino = os.path.join(destino, f"{pdb_id}.pdb")

    try:
        print(f"Descargando estructura PDB '{pdb_id}' desde {base_url}...")
        response = requests.get(base_url, timeout=10)
        response.raise_for_status()  # Lanza una excepción si la descarga falla

        with open(archivo_destino, 'w') as f:
            f.write(response.text)

        print(f"Archivo PDB '{pdb_id}' descargado exitosamente en '{archivo_destino}'.")
        return archivo_destino
    except requests.RequestException as e:
        print(f"Error al descargar la estructura PDB '{pdb_id}': {e}")
        return None

def DescargarArchivo():
    """
    Descarga un archivo PDB y lo filtra para eliminar líneas no deseadas ('HETATM' y 'CONECT').

    Utiliza variables globales 'pdb_id' y 'filename' para almacenar información sobre el archivo.
    """
    global filename
    global pdb_id

    archivo_destino = os.getcwd()
    filename = descargar_pdb_directo(pdb_id, destino=archivo_destino)

    if not filename:
        print("Could not download the file. Please check the PDB ID and your connection.")
        return

    filtered_lines = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith("HETATM") and not line.startswith("CONECT"):
                    filtered_lines.append(line)

        with open(filename, 'w') as f:
            f.writelines(filtered_lines)

        print(f"PDB file '{pdb_id}' processed and saved as '{filename}'.")
    except Exception as e:
        print(f"Error processing the file: {e}")
    
##Centrar la ventana root
def center_window(root, width, height):
    # Obtener las dimensiones de la pantalla
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()
    
    # Calcular las coordenadas para centrar la ventana en el lado izquierdo
    x = 0
    y = (screen_height - height) // 2  # Centrar verticalmente
    
    # Establecer la geometría de la ventana
    root.geometry(f"{width}x{height}+{x}+{y}")


def colornodo(event):
    global col_n
    col_n = colores_n.get()
    
def colorarista(event):
    global col_a
    col_a = colores_a.get()

def coloretiqueta(event):
    global col_e
    col_e = colores_l.get()

def Distancia():
    global distancia
    global dis
    dis = distancia.get()
    dis = float(dis)

def Distancia2():
    global distancia2
    global dis2
    dis2 = distanciaP2.get()
    dis2 = float(dis2)

def IDPDB():
    global pdb_id
    pdb_id= id_pdb.get()

def actualizar_estado_color():
    global op_lab
    op_lab = opcion_la.get()
    if op_lab == "Yes":
        colores_l.config(state="readonly")
    else:
        colores_l.set("")  # Reiniciar la selección
        colores_l.config(state="disabled")


def ConstruirRED():
    global distancia
    global filename
    global col_n
    global col_a
    global col_e
    global op_lab
    if filename is None:
        print("First, select a file.")
        return
    warnings.filterwarnings("ignore", category=UserWarning, module="Bio.PDB")
    parser=PDBParser()
    p = PDBParser()
    structure = p.get_structure("X",filename)
    j = 0

    resi={}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    resi[j]=residue
                    j=j+1
    n=len(resi)            
    posiciones = {}
    for nodo in resi:
        posiciones[nodo]=tuple(resi[nodo]["CA"].get_coord())
    G=nx.Graph()

    for nodo in posiciones:
        G.add_node(nodo,pos=posiciones[nodo])
   #Color nodos 
    if col_n == "White":
        rgb_n = (1.0, 1.0, 1.0)
    elif col_n == "Black":
        rgb_n = (0.0, 0.0, 0.0)
    elif col_n == "Red":
        rgb_n = (1.0, 0.0, 0.0)
    elif col_n == "Green":
        rgb_n = (0.0, 1.0, 0.0)
    elif col_n == "Blue":
        rgb_n = (0.0, 0.0, 1.0)
    elif col_n == "Yellow":
        rgb_n = (1.0, 1.0, 0.0)
    elif col_n == "Magenta":
        rgb_n = (1.0, 0.0, 1.0)
    elif col_n == "Cyan":
        rgb_n = (0.0, 1.0, 1.0)
    elif col_n == "Gray":
        rgb_n = (0.3, 0.3, 0.3)
    elif col_n == "Orange":
        rgb_n = (1.0, 0.65, 0.0)
    elif col_n == "Blue":
        rgb_n = (0.0, 0.0, 1.0)
    #Color aristas    
    if col_a == "White":
        rgb_a = (1.0, 1.0, 1.0)
    elif col_a == "Black":
        rgb_a = (0.0, 0.0, 0.0)
    elif col_a == "Red":
        rgb_a = (1.0, 0.0, 0.0)
    elif col_a == "Green":
        rgb_a = (0.0, 1.0, 0.0)
    elif col_a == "Blue":
        rgb_a = (0.0, 0.0, 1.0)
    elif col_a == "Yellow":
        rgb_a = (1.0, 1.0, 0.0)
    elif col_a == "Magenta":
        rgb_a = (1.0, 0.0, 1.0)
    elif col_a == "Cyan":
        rgb_a = (0.0, 1.0, 1.0)
    elif col_a == "Gray":
        rgb_a = (0.3, 0.3, 0.3)
    elif col_a == "Orange":
        rgb_a = (1.0, 0.65, 0.0)
    elif col_a == "Blue":
        rgb_a = (0.0, 0.0, 1.0)
    #Colores etiqueta
    if col_e == ("Black"):
        rgb_e = (0.0,0.0,0.0)
    elif col_e == ("White"):
        rgb_e = (1.0, 1.0, 1.0)
                


    for i in range (n-1):
        for j in range (i+1,n):
            if resi[i]["CA"]-resi[j]["CA"]<=dis:
                G.add_edges_from([(i,j)])
    tr = {}
    num = {}
    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        numero_residuo = residuo.id[1]
        tr[numero_res]=tipo_residuo
        num[numero_res]=numero_residuo
        
    if op_lab =="Yes":
        for node, (x, y, z) in posiciones.items():
            mlab.points3d(x, y, z, color=(rgb_n), scale_factor=1.0)
            etiqueta = f"{tr[node]}{num[node]}"
            mlab.text3d(x, y, z, etiqueta, scale=(0.4, 0.4, 0.4), color=rgb_e)
    else:
        for node, (x, y, z) in posiciones.items():
            mlab.points3d(x, y, z, color=(rgb_n), scale_factor=1.0)
            
        
    for edge in G.edges():
        i, j = edge
        x, y, z = zip(*[posiciones[i], posiciones[j]])
        mlab.plot3d(x, y, z, tube_radius=0.1, color=rgb_a)
    


    # Ajusta la vista 3D
    mlab.view(azimuth=0, elevation=90, distance=5)

    # Muestra la visualización interactiva
    mlab.show()

def extraer_nombre_proteina(pdb_path):
    nombre_proteina = None

    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            # Buscar la línea con el prefijo TITLE
            if line.startswith("TITLE"):
                if nombre_proteina is None:
                    nombre_proteina = line[10:].strip()
                else:
                    # Continuar el nombre si hay múltiples líneas TITLE
                    nombre_proteina += " " + line[10:].strip()

            # Terminar después del encabezado si ya se ha leído el nombre
            if line.startswith("ATOM"):
                break

    return nombre_proteina

def propiedad_Funcion(event):
    global prop
    prop = propiedades_lista.get()

def lambda_dos(G):
    centralidad = nx.eigenvector_centrality(G, max_iter=1000)
    A = nx.to_numpy_array(G)
    x=np.array(list(centralidad.values()))
    Ax=A @ x
    
    lambda_value= Ax/x
    lambda_valores = np.linalg.eigvals(A)
    orden_lambda = np.sort(lambda_valores)[::-1]
    maximos_lamda= orden_lambda[:2]
    m1 = maximos_lamda[0]
    m2 = maximos_lamda[1]
    tag1=np.where(lambda_valores == m1)[0]
    tag2=np.where(lambda_valores == m2)[0]
    
    I= np.eye(A.shape[0])
    M= A- m2* I
    _, v= np.linalg.eig(A)
    vector_propio2=v[:,np.isclose(_,m2)]
    
    vector_propio2_real = np.real(vector_propio2).flatten()
    
    # Crear un diccionario con los valores reales del vector propio
    lambda_2_dict = {nodo: vector_propio2_real[i] for i, nodo in enumerate(G.nodes())}
    
    return lambda_2_dict

def Propiedad_Centralidad():
    global filename
    global prop
    global dis2
    if filename is None:
        print("First, select a file.")
        return
    warnings.filterwarnings("ignore", category=UserWarning, module="Bio.PDB")
    parser=PDBParser()
    p = PDBParser()
    structure = p.get_structure("X",filename)
    j = 0

    resi={}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    resi[j]=residue
                    j=j+1
    n=len(resi)            
    posiciones = {}
    for nodo in resi:
        posiciones[nodo]=tuple(resi[nodo]["CA"].get_coord())
        G=nx.Graph()

    for nodo in posiciones:
        G.add_node(nodo,pos=posiciones[nodo])
              


    for i in range (n-1):
        for j in range (i+1,n):
            if resi[i]["CA"]-resi[j]["CA"]<=dis2:
                G.add_edges_from([(i,j)])   
    
    if prop == "Degree":
        # Calcular el grado de cada nodo
        prop_cent = dict(G.degree())
        pdc = "Degree" 
       
    elif prop =="Betweenness":
        prop_cent = dict(nx.betweenness_centrality(G))
        pdc = "Betweenness"
                        
    elif prop == "Closeness":
        prop_cent = dict(nx.closeness_centrality(G))
        pdc = "Closeness"
    
                        
    elif prop == "Clustering Coefficient":
        prop_cent = dict(nx.clustering(G))
        pdc = "Clustering Coefficient"
        
                       
    elif prop == "Eigenvector":
        prop_cent = dict(nx.eigenvector_centrality(G, max_iter=1000))
        pdc = "Eigenvector"
        
    elif prop == "Eigenvector Lamda2":
        prop_cent = dict(lambda_dos(G))
        pdc = "Eigenvector Lamda2"
        
    elif prop == "Excentricidad":
        prop_cent = dict(nx.eccentricity(G))
        pdc = "Eccentricity"
        
    elif prop ==  "Stress":
        prop_cent = dict(Stress_centrality(G))
        pdc =  "Stress"
        
                       

    
                                    
    # Normalizar los valores 
    min_prop_cent = min(prop_cent.values())
    max_prop_cent = max(prop_cent.values())
    norm_prop_cent = {nodo: (valor - min_prop_cent) / (max_prop_cent - min_prop_cent) for nodo, valor in prop_cent.items()}

    # Definir colores según los rangos especificados
    colores = np.zeros((len(G), 3))  # Inicializar matriz de colores
    for i, (nodo, prop_cent_normalizado) in enumerate (norm_prop_cent.items()):
        if prop_cent_normalizado == 1:
            colores[i] = [1.0, 0.0, 0.0]  # Rojo intenso
        elif 0.8 <= prop_cent_normalizado < 1:
            colores[i] = [1.0, 0.325, 0.0]  # Rojo anaranjado
        elif 0.6 <= prop_cent_normalizado < 0.8:
            colores[i] = [1.0, 0.647, 0.0]  # Naranja
        elif 0.4 <= prop_cent_normalizado < 0.6:
            colores[i] = [1.0, 0.745, 0.471]  # Naranja suave
        elif 0.2 <= prop_cent_normalizado < 0.4:
            colores[i] = [1.0, 0.843, 0.647]  # Color carnita o durazno
        elif 0.1 <= prop_cent_normalizado < 0.2:
            colores[i] = [1.0, 0.941, 0.824]  # Durazno más claro
        elif 0 <= prop_cent_normalizado < 0.1:
            colores[i] = [0.5, 0.5, 0.5]  # Gris


    tr = {}
    num = {}
    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        numero_residuo = residuo.id[1]
        tr[numero_res]=tipo_residuo
        num[numero_res]=numero_residuo

            

    dfcentralidad = pd.DataFrame()
    propiedades = ["Node", "Amino acid", pdc]
    datos = []
            
    for i in range(len(tr)):
        datos.append({
        "Node": i,
        "Amino acid":  f"{tr[i]} {num[i]}"
        }) 
           
    dfcentralidad = pd.DataFrame(datos, columns=propiedades)
    prop_cent_series = pd.Series(prop_cent)

    # Agregar los valores a una nueva columna en dfDatos
    dfcentralidad[pdc] = prop_cent_series

    pd.set_option('display.max_rows', None)
    print(dfcentralidad)
    
    ######################################################################################
    # Agrega nodos con colores
    for i, nodo in enumerate(G.nodes):
        x, y, z = posiciones[nodo]
        color = tuple(colores[i])  # Convertir la matriz a tupla
        etiqueta = f"{tr[nodo]}{num[nodo]}"
        mlab.text3d(x, y, z, etiqueta, scale=(0.5, 0.5, 0.5), color=(1,1,1),)
        mlab.points3d(x, y, z, color=color, scale_factor=1.0)

    mlab.text(0.5, 0.9, f"{pdc}", width=0.2, color=(1, 1, 1), line_width=1.5)
    # Agrega bordes
    for edge in G.edges:
        i, j = edge
        x, y, z = zip(posiciones[i], posiciones[j])
        mlab.plot3d(x, y, z, tube_radius=0.05, color=(0.5, 0.5, 0.5))

    # Muestra la visualización
    mlab.show()

#dfcentralidad_excel = None

def MostrarPropiedadesS():
    global dfcetralidad_excel
    global filename
    global prop
    global distancia2
    dfcentralidad_excel = pd.DataFrame() 
    if filename is None:
        print("First, select a file.")
        return

    parser=PDBParser()
    p = PDBParser()
    structure = p.get_structure("X",filename)
    j = 0

    resi={}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    resi[j]=residue
                    j=j+1
    n=len(resi)            
    posiciones = {}
    for nodo in resi:
        posiciones[nodo]=tuple(resi[nodo]["CA"].get_coord())
        G=nx.Graph()

    for nodo in posiciones:
        G.add_node(nodo,pos=posiciones[nodo])
              


    for i in range (n-1):
        for j in range (i+1,n):
            if resi[i]["CA"]-resi[j]["CA"]<=dis2:
                G.add_edges_from([(i,j)])   
    
    if prop == "Degree":
        # Calcular el grado de cada nodo
        prop_cent = dict(G.degree())
        pdc = "Degree" 
       
    elif prop =="Betweenness":
        prop_cent = dict(nx.betweenness_centrality(G))
        pdc = "Betweenness"
                        
    elif prop == "Closeness":
        prop_cent = dict(nx.closeness_centrality(G))
        pdc = "Closeness"
    
                        
    elif prop == "Clustering Coefficient":
        prop_cent = dict(nx.clustering(G))
        pdc = "Clustering Coefficient"
        
                       
    elif prop == "Eigenvector":
        prop_cent = dict(nx.eigenvector_centrality(G, max_iter=1000))
        pdc = "Eigenvector"
        
    elif prop == "Eigenvector Lamda2":
        prop_cent = dict(lambda_dos(G))
        pdc = "Eigenvector Lamda2"
        
    elif prop ==  "Eccentricity":
        prop_cent = dict(nx.eccentricity(G))
        pdc =  "Eccentricity"
        
    elif prop == "Stress":
        prop_cent = dict(Stress_centrality(G))
        pdc = "Stress"
                        
    tr = {}
    num = {}
    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        numero_residuo = residuo.id[1]
        tr[numero_res]=tipo_residuo
        num[numero_res]=numero_residuo

            

    dfcentralidad = pd.DataFrame()
    propiedades = ["Node", "Amino acid", pdc]
    datos = []
            
    for i in range(len(tr)):
        datos.append({
        "Node": i,
        "Amino acid":  f"{tr[i]} {num[i]}"
        }) 
           
    dfcentralidad = pd.DataFrame(datos, columns=propiedades)
    if isinstance(prop_cent, dict):
        prop_cent_series = pd.Series(prop_cent)
    else:
        # Si prop_cent es un array o lista, extrae los valores escalares
        prop_cent_series = pd.Series([x[0] if isinstance(x, (list, np.ndarray)) else x for x in prop_cent.values()])

    

    # Agregar los valores a una nueva columna en dfDatos
    dfcentralidad[pdc] = prop_cent_series

    plt.plot(dfcentralidad["Amino acid"], dfcentralidad[pdc], marker='o')
    plt.xlabel("Amino acid")
    plt.ylabel(pdc)
    plt.title(f"Grafico de {pdc}")
    plt.show

    if pdc != "Eigenvector Lamda2":
        num_nodos = G.number_of_nodes()
        num_aristas = G.number_of_edges()
        densidad = nx.density(G)
        nombre = extraer_nombre_proteina(filename)
    
        # Crear un DataFrame con estas propiedades
        propiedades_red = {
            "Property": ["Protein name", "Number of nodes", "Number of edges", "Density"],
            "Value": [nombre, num_nodos, num_aristas, densidad],
        }
        df_resumen = pd.DataFrame(propiedades_red)
        df_resumen_str = df_resumen.to_string(index = False)
        
        pd.set_option('display.max_rows', None)
        dfcentralidad_excel=dfcentralidad.copy()
        #print(dfcentralidad)
        
        # Convertir el DataFrame a una cadena de texto
        dataframe_str = dfcentralidad.to_string(index=False)
    
        # Limpiar el contenido actual del widget Text
        dataframe_text.delete(1.0, tk.END)
        
        dataframe_text.insert (tk.END, df_resumen_str)
        
        dataframe_text.insert(tk.END, "\n\n")
        
        promedio = dfcentralidad[pdc].mean()
        
        
        valor_maximo = dfcentralidad[pdc].max()
        filas_maximas = dfcentralidad[dfcentralidad[pdc] == valor_maximo]
        
        valor_minimo = dfcentralidad[pdc].min()
        filas_minimas = dfcentralidad[dfcentralidad[pdc] == valor_minimo]
        
        # Insertar el DataFrame en el widget Text
        dataframe_text.insert(tk.END, f"Maximum value:\n {filas_maximas}\n\n")
        dataframe_text.insert(tk.END, f"Minimum value:\n {filas_minimas}\n\n")
        dataframe_text.insert(tk.END, f"The average of {pdc} is {promedio}\n\n")
        dataframe_text.insert(tk.END, dataframe_str)
    else:
        pd.set_option('display.max_rows', None)
        dfcentralidad_excel=dfcentralidad.copy()
        #print(dfcentralidad)
        
        # Convertir el DataFrame a una cadena de texto
        dataframe_str = dfcentralidad.to_string(index=False)
    
        # Limpiar el contenido actual del widget Text
        dataframe_text.delete(1.0, tk.END)
    
        # Insertar el DataFrame en el widget Text
        dataframe_text.insert(tk.END, dataframe_str)
    
    
    # Crear una ventana de Tkinter para mostrar el gráfico
    grafico_ventana = tk.Toplevel()
    grafico_ventana.title(f"{pdc} plot")

    # Crear un objeto de figura
    fig = Figure(figsize=(30, 8))
    ax = fig.add_subplot(111)
    ax.plot(dfcentralidad["Node"], dfcentralidad[pdc])
    ax.set_xlabel("Nodes")    
    ax.set_ylabel(pdc)
    ax.set_title(f"{pdc} plot")
    ax.tick_params(axis='x', rotation=90)
    ax.axhline(0, color='black', linestyle='--', linewidth=1, label='División y=0')
    ax.grid()
    
    # Mostrar la figura en la ventana de Tkinter
    canvas = FigureCanvasTkAgg(fig, master=grafico_ventana)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    def guardar_grafico():
        # Abrir un cuadro de diálogo para guardar archivo
        archivo = tk.filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
        )
        if archivo:
            fig.savefig(archivo)
            print(f"Graph saved at: {archivo}")

    # Botón para guardar el gráfico
    boton_guardar = tk.Button(grafico_ventana, text="Save plot", command=guardar_grafico)
    boton_guardar.pack()


    # Agregar un botón para cerrar la ventana
    boton_cerrar = tk.Button(grafico_ventana, text="Close", command=grafico_ventana.destroy)
    boton_cerrar.pack()

def M_T_V():
    global dfcetralidad_excel
    global filename
    global prop
    global dis2
    if filename is None:
        print("First, select a file.")
        return
    

    parser=PDBParser()
    p = PDBParser()
    structure = p.get_structure("X",filename)
    j = 0

    resi={}
    cadena={}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    resi[j]=residue
                    cadena[j]=chain
                    j=j+1
    n=len(resi)            
    posiciones = {}
    for nodo in resi:
        posiciones[nodo]=tuple(resi[nodo]["CA"].get_coord())
        G=nx.Graph()

    for nodo in posiciones:
        G.add_node(nodo,pos=posiciones[nodo])
              


    for i in range (n-1):
        for j in range (i+1,n):
            if resi[i]["CA"]-resi[j]["CA"]<=dis2:
                G.add_edges_from([(i,j)])   
    tr = {}
    num = {}
   
    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        numero_residuo = residuo.id[1]
        tr[numero_res]=tipo_residuo
        num[numero_res]=numero_residuo
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Propiedades de centralidad~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    grado=nx.degree(G)
    intermediacion=nx.betweenness_centrality(G)
    cercania=nx.closeness_centrality(G)
    agrupamiento=nx.clustering(G)
    vector_propio = nx.eigenvector_centrality(G, max_iter=1000)
    vector_propio_lambda2 = lambda_dos(G)
    excentricidad = nx.eccentricity(G)
    estres_v = Stress_centrality(G)
    
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tipo de aminoacido~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aminoacidos_no_polares = {'ALA', 'VAL', 'ILE', 'LEU', 'PHE', 'MET', 'PRO', 'TRP', 'GLY'}
    aminoacidos_polar = {'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'}
    aminoacidos_positivos = {'ARG', 'LYS', 'HIS'}
    aminoacidos_negativos = {'ASP', 'GLU'}

    aminoacidos_alifaticos = {'ALA', 'VAL', 'ILE', 'LEU'}
    aminoacidos_aromaticos = {'PHE', 'TYR', 'TRP'}
    aminoacidos_con_azufre = {'CYS', 'MET'}
    aminoacidos_con_grupos_hidroxilo = {'SER', 'THR'}
    aminoacidos_basicos = {'ARG', 'LYS', 'HIS'}
    aminoacidos_acidos = {'ASP', 'GLU'}

    
    
    n_resi={}

    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        n_resi[numero_res]=tipo_residuo

    tipo_de_aminoacido ={}

    for r in range (len(n_resi)):
        if n_resi[r] in aminoacidos_no_polares:
            tipo_de_aminoacido[r] = "Non-polar"   
        elif n_resi[r] in aminoacidos_polar:
            tipo_de_aminoacido[r] = "Polar"
        elif n_resi[r] in aminoacidos_positivos:
            tipo_de_aminoacido[r] = "Positive"
        elif n_resi[r] in aminoacidos_negativos:
            tipo_de_aminoacido[r] = "Negative"
        else:
            tipo_de_aminoacido[r] =  "Other"
            
            
    tipo_de_aminoacido_G ={}

    for r_G in range (len(n_resi)):
        if n_resi[r_G] in aminoacidos_alifaticos:
            tipo_de_aminoacido_G[r] = "Non-polar"   
        elif n_resi[r_G] in aminoacidos_aromaticos:
            tipo_de_aminoacido_G[r] = "Aromatic"
        elif n_resi[r_G] in aminoacidos_con_azufre:
            tipo_de_aminoacido_G[r] = "With Sulfur"
        elif n_resi[r_G] in aminoacidos_con_grupos_hidroxilo:
            tipo_de_aminoacido_G[r_G] = "With Hydroxyl group"
        elif n_resi[r_G] in aminoacidos_basicos:
            tipo_de_aminoacido_G[r_G] = "Basic"
        elif n_resi[r_G] in aminoacidos_acidos:
            tipo_de_aminoacido_G[r_G] = "Acidic"
        else:
            tipo_de_aminoacido_G[r_G] =  "Other"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            
    dfcen_col_dis=pd.DataFrame()
    propiedades = ["Node", "Amino acid", "Chain", "Amino acid type", "Degree", "Betweenness", "Closeness", "Clustering Coefficient", "Eigenvector", "Eigenvector Lamda2", "Eccentricity", "Stress"]
    datos = []
    #Si queremos montrar el numero de aminoacidos en la lista sumar el desface en "1.Nodo": i
    for i in range(len(nx.clustering(G))):
        datos.append({
            "Node": i+1,
            "Amino acid":  f"{tr[i]} {num[i]}",
            "Chain": str(cadena[i]) ,
            "Amino acid type": tipo_de_aminoacido[i],
            "Degree": grado[i],
            "Betweenness": intermediacion[i],
            "Closeness": cercania[i],
            "Clustering Coefficient": agrupamiento[i],
            "Eigenvector": vector_propio[i],
            "Eigenvector Lamda2": vector_propio_lambda2[i],
            "Eccentricity": excentricidad[i],
            "Stress": estres_v[i]
        })

    dfcen_col_dis = pd.DataFrame(datos, columns=propiedades)
    print(len(G.edges))
    #dfcen_col_dis.to_excel("Prop_mich0-7.xlsx", index=False)
    #print (dfcen_col_dis)
    dfcentralidad_excel=dfcen_col_dis
    dataframe_str = dfcen_col_dis.to_string(index=False)
    
    num_nodos = G.number_of_nodes()
    num_aristas = G.number_of_edges()
    densidad = nx.density(G)
    nombre = extraer_nombre_proteina(filename)

    # Crear un DataFrame con estas propiedades
    propiedades_red = {
        "Property": ["Protein name", "Number of nodes", "Number of edges", "Density"],
        "Value": [nombre, num_nodos, num_aristas, densidad],
    }
    df_resumen = pd.DataFrame(propiedades_red)
    df_resumen_str = df_resumen.to_string(index = False)
    
    # Limpiar el contenido actual del widget Text
    dataframe_text.delete(1.0, tk.END)

    
    dataframe_text.insert (tk.END, df_resumen_str)
    dataframe_text.insert(tk.END, "\n\n")


    # Insertar el DataFrame en el widget Text
    dataframe_text.insert(tk.END, dataframe_str)

# Función para crear la ventana emergente con los Checkbox
def seleccionar_propiedades():
    # Crear una ventana emergente
    ventana_seleccion = tk.Toplevel()
    ventana_seleccion.title("Select properties")
    
    # Lista de propiedades de centralidad
    propiedades = ["Degree", "Betweenness", "Closeness", "Clustering Coefficient", "Eigenvector", "Eigenvector Lamda2", "Eccentricity", "Stress"]
    
    # Diccionario para almacenar los valores de los Checkbox
    seleccion = {prop: tk.BooleanVar() for prop in propiedades}
    
    # Crear un Checkbox para cada propiedad
    for i, prop in enumerate(propiedades):
        tk.Checkbutton(ventana_seleccion, text=prop, variable=seleccion[prop]).grid(row=i, column=0, sticky="w")
    
    # Botón para confirmar la selección
    tk.Button(ventana_seleccion, text="Accept", command=ventana_seleccion.destroy).grid(row=len(propiedades), column=0, pady=10)
    
    # Esperar a que la ventana se cierre
    ventana_seleccion.wait_window()
    
    # Retornar las propiedades seleccionadas
    return [prop for prop, var in seleccion.items() if var.get()]

# Función modificada para guardar los datos
def Guardar():
    global dfcentralidad_excel
    global filename
    global prop
    global dis2
    
    if filename is None:
        print("First, select a file.")
        return
    
    # Abrir la ventana para seleccionar propiedades
    propiedades_seleccionadas = seleccionar_propiedades()
    
    if not propiedades_seleccionadas:
        messagebox.showinfo("Cancelled", "No properties selected.")
        return
    
    # Procesar el archivo PDB
    parser = PDBParser()
    p = PDBParser()
    structure = p.get_structure("X", filename)
    j = 0

    resi = {}
    cadena = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    resi[j] = residue
                    cadena[j] = chain
                    j = j + 1
    n = len(resi)
    posiciones = {}
    for nodo in resi:
        posiciones[nodo] = tuple(resi[nodo]["CA"].get_coord())
        G = nx.Graph()

    for nodo in posiciones:
        G.add_node(nodo, pos=posiciones[nodo])

    for i in range(n - 1):
        for j in range(i + 1, n):
            if resi[i]["CA"] - resi[j]["CA"] <= dis2:
                G.add_edges_from([(i, j)])
    
    tr = {}
    num = {}
    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        numero_residuo = residuo.id[1]
        tr[numero_res] = tipo_residuo
        num[numero_res] = numero_residuo

    # Calcular solo las propiedades seleccionadas
    resultados = {}
    if "Degree" in propiedades_seleccionadas:
        resultados["Degree"] = nx.degree(G)
    if "Betweenness" in propiedades_seleccionadas:
        resultados["Betweenness"] = nx.betweenness_centrality(G)
    if "Closeness" in propiedades_seleccionadas:
        resultados["Closeness"] = nx.closeness_centrality(G)
    if "Clustering Coefficient" in propiedades_seleccionadas:
        resultados["Clustering Coefficient"] = nx.clustering(G)
    if "Eigenvector" in propiedades_seleccionadas:
        resultados["Eigenvector"] = nx.eigenvector_centrality(G, max_iter=1000)
    if "Eigenvector Lamda2" in propiedades_seleccionadas:
        resultados["Eigenvector Lamda2"] = lambda_dos(G)
    if "Eccentricity" in propiedades_seleccionadas:
        resultados["Eccentricity"] = nx.eccentricity(G)
    if "Stress" in propiedades_seleccionadas:
        resultados["Stress"] = Stress_centrality(G)
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tipo de aminoacido~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    aminoacidos_no_polares = {'ALA', 'VAL', 'ILE', 'LEU', 'PHE', 'MET', 'PRO', 'TRP', 'GLY'}
    aminoacidos_polar = {'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'}
    aminoacidos_positivos = {'ARG', 'LYS', 'HIS'}
    aminoacidos_negativos = {'ASP', 'GLU'}

    aminoacidos_alifaticos = {'ALA', 'VAL', 'ILE', 'LEU'}
    aminoacidos_aromaticos = {'PHE', 'TYR', 'TRP'}
    aminoacidos_con_azufre = {'CYS', 'MET'}
    aminoacidos_con_grupos_hidroxilo = {'SER', 'THR'}
    aminoacidos_basicos = {'ARG', 'LYS', 'HIS'}
    aminoacidos_acidos = {'ASP', 'GLU'}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Resumen~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    num_nodos = G.number_of_nodes()
    num_aristas = G.number_of_edges()
    densidad = nx.density(G)
    nombre = extraer_nombre_proteina(filename)

    # Crear un DataFrame con estas propiedades
    propiedades_red = {
        "Property": ["Protein name", "Number of nodes", "Number of edges", "Density"],
        "Value": [nombre, num_nodos, num_aristas, densidad],
    }
    df_resumen = pd.DataFrame(propiedades_red)

    

    
    
    n_resi={}

    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        n_resi[numero_res]=tipo_residuo

    tipo_de_aminoacido ={}

    for r in range (len(n_resi)):
        if n_resi[r] in aminoacidos_no_polares:
            tipo_de_aminoacido[r] = "Non-polar"   
        elif n_resi[r] in aminoacidos_polar:
            tipo_de_aminoacido[r] = "Polar"
        elif n_resi[r] in aminoacidos_positivos:
            tipo_de_aminoacido[r] = "Positive"
        elif n_resi[r] in aminoacidos_negativos:
            tipo_de_aminoacido[r] = "Negative"
        else:
            tipo_de_aminoacido[r] =  "Other"
            
            
    tipo_de_aminoacido_G ={}

    for r_G in range (len(n_resi)):
        if n_resi[r_G] in aminoacidos_alifaticos:
            tipo_de_aminoacido_G[r] = "Non-polar"   
        elif n_resi[r_G] in aminoacidos_aromaticos:
            tipo_de_aminoacido_G[r] = "Aromatic"
        elif n_resi[r_G] in aminoacidos_con_azufre:
            tipo_de_aminoacido_G[r] = "With Sulfur"
        elif n_resi[r_G] in aminoacidos_con_grupos_hidroxilo:
            tipo_de_aminoacido_G[r_G] = "With Hydroxyl group"
        elif n_resi[r_G] in aminoacidos_basicos:
            tipo_de_aminoacido_G[r_G] = "Basic"
        elif n_resi[r_G] in aminoacidos_acidos:
            tipo_de_aminoacido_G[r_G] = "Acidic"
        else:
            tipo_de_aminoacido_G[r_G] =  "Other"
    

    # Crear el DataFrame con las propiedades seleccionadas
    datos = []
    for i in range(len(nx.clustering(G))):
        fila = {
            "Node": i + 1,
            "Amino acid": f"{tr[i]} {num[i]}",
            "Chain": str(cadena[i]),
            "Amino acid type": tipo_de_aminoacido[i]
        }
        for prop in propiedades_seleccionadas:
            fila[prop] = resultados[prop][i]
        datos.append(fila)
    
    dfcen_col_dis = pd.DataFrame(datos, columns=["Node", "Amino acid", "Chain", "Amino acid type"] + propiedades_seleccionadas)

    # Preguntar por el formato deseado
    formato = simpledialog.askstring(
        "File format",
        "What format would you like to save the file in? (xlsx, csv, txt, json)"
    )

    if not formato:
        messagebox.showinfo("Cancelled", "No file format selected.")
        return

    # Abrir cuadro de diálogo para guardar el archivo
    file_path = filedialog.asksaveasfilename(
        defaultextension=f".{formato}",
        filetypes=[
            ("Excel files", "*.xlsx"),
            ("CSV files", "*.csv"),
            ("Text files", "*.txt"),
            ("JSON files", "*.json"),
            ("All files", "*.*")
        ],
        title="Save file as"
    )

    if not file_path:
        messagebox.showinfo("Cancelled", "No file selected to save.")
        return

    # Guardar el archivo según el formato elegido
    try:
        if formato.lower() == "xlsx":
            dfcen_col_dis.to_excel(file_path, index=False)
        elif formato.lower() == "csv":
            dfcen_col_dis.to_csv(file_path, index=False)
        elif formato.lower() == "txt":
            dfcen_col_dis.to_csv(file_path, index=False, sep='\t')
        elif formato.lower() == "json":
            dfcen_col_dis.to_json(file_path, orient='records', indent=4)
        else:
            messagebox.showerror("Error", f"Unsupported format: {formato}")
            return

        messagebox.showinfo("Success", f"The file has been successfully saved as {formato}.")
    except Exception as e:
        messagebox.showerror("Error", f"Could not save the file: {e}")

def DistanciasRed():
    global filename
    if filename is None:
        print("First, select a file.")
        return
    parser=PDBParser()
    p = PDBParser()
    structure = p.get_structure("X",filename)
    j = 0

    resi={}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    resi[j]=residue
                    j=j+1
    n=len(resi)            
    # Creamos un diccionario para almacenar las posiciones de los carbonos alfa
    posiciones = {}

    for nodo in resi:
        posiciones[nodo]=tuple(resi[nodo]["CA"].get_coord())


    # Creamos una red "G" y añadimos los nodos de nuestra lista resi
    G=nx.Graph()
    for nodo in posiciones:
        G.add_node(nodo,pos=posiciones[nodo])

    # Creamos los bordes de la red y asignamos color
    for i in range (n-1):
        for j in range (i+1,n):
            distancia=resi[i]["CA"]-resi[j]["CA"]
            if distancia>=3.8 and distancia <=5:
                G.add_edges_from([(i,j, {'color': 'blue'})])
            elif distancia>5 and distancia <=6:
                G.add_edges_from([(i,j, {'color': 'yellow'})])
            elif distancia>6 and distancia <=7:
                G.add_edges_from([(i,j, {'color': 'red'})])

    # Generar grafico 3D
    tr = {}
    num = {}
    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        numero_residuo = residuo.id[1]
        tr[numero_res]=tipo_residuo
        num[numero_res]=numero_residuo

    for edge in G.edges(data=True):
        i, j, color_data = edge
        x, y, z = zip(*[posiciones[i], posiciones[j]])
        color = color_data['color']  # Obtiene el color desde los datos de la arista
        if color == 'blue':
            color = (0, 0, 1)  
        elif color == 'yellow':
            color = (1, 1, 0)  
        elif color == 'red':
            color = (1, 0, 0)  
        mlab.plot3d(x, y, z, tube_radius=0.1, color=color)

    for node, (x, y, z) in posiciones.items():
        mlab.points3d(x, y, z, color=((0, 1, 0)), scale_factor=1.0)
        etiqueta = f"{tr[node]}{num[node]}"
        mlab.text3d(x, y, z, etiqueta, scale=(0.4, 0.4, 0.4), color=(0,0, 0),)
    # Ajusta la vista 3D
    mlab.view(azimuth=0, elevation=90, distance=5)

    # Muestra la visualización interactiva
    mlab.show()

def Rangos_Dis1_1(event):
    global dis1_1
    dis1_1 = float(distancia1_1.get())


def Rangos_Dis1_2(event):    
    global dis1_2
    dis1_2 = float(distancia1_2.get())

    
def Rangos_Dis2_1(event):
    global dis2_1
    dis2_1 = float(distancia2_1.get())
    
def Rangos_Dis2_2(event):
    global dis2_2
    dis2_2 = float(distancia2_2.get())
    
def Rangos_Dis3_1(event):
    global dis3_1
    dis3_1 = float(distancia3_1.get())
    
def Rangos_Dis3_2(event):
    global dis3_2
    dis3_2 = float(distancia3_2.get())
    
def Colores_Rango1(event=None):
    global colorR1
    
    colorR1 = ColorRango1.get()
    print(colorR1)
   
    
def Colores_Rango2(event=None):
    global colorR2
    
    colorR2 = ColorRango2.get()
    print(colorR2)
    
def Colores_Rango3(event=None):
    global colorR3

    colorR3 = ColorRango3.get()
    print(colorR3)

def Colores_Rango_Nodo(event=None):
    global colorRNodo

    colorRNodo = ColorNodo.get()
    print(colorRNodo)


def DistanciasRedPersonalizada():
    global filename
    global dis1_1
    global dis1_2
    global dis2_1
    global dis2_2
    global dis3_1
    global dis3_2
    global colorR1
    global colorR2
    global colorR3
    global colorRNodo
    if filename is None:
        print("First, select a file.")
        return
    parser=PDBParser()
    p = PDBParser()
    structure = p.get_structure("X",filename)
    j = 0

    resi={}
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    resi[j]=residue
                    j=j+1
    n=len(resi)            
    # Creamos un diccionario para almacenar las posiciones de los carbonos alfa
    posiciones = {}

    for nodo in resi:
        posiciones[nodo]=tuple(resi[nodo]["CA"].get_coord())


    # Creamos una red "G" y añadimos los nodos de nuestra lista resi
    G=nx.Graph()
    for nodo in posiciones:
        G.add_node(nodo,pos=posiciones[nodo])

    # Creamos los bordes de la red y asignamos color
    for i in range (n-1):
        for j in range (i+1,n):
            distancia=resi[i]["CA"]-resi[j]["CA"]
            if distancia>=dis1_1 and distancia <=dis1_2:
                G.add_edges_from([(i,j, {'color': colorR1})])
            elif distancia>dis2_1 and distancia <=dis2_2:
                G.add_edges_from([(i,j, {'color': colorR2})])
            elif distancia>dis3_1 and distancia <=dis3_2:
                G.add_edges_from([(i,j, {'color': colorR3})])

    # Generar grafico 3D
    tr = {}
    num = {}
    for numero_res, residuo in resi.items():
        tipo_residuo = residuo.resname
        numero_residuo = residuo.id[1]
        tr[numero_res]=tipo_residuo
        num[numero_res]=numero_residuo

    for edge in G.edges(data=True):
        i, j, color_data = edge
        x, y, z = zip(*[posiciones[i], posiciones[j]])
        color = color_data['color']  # Obtiene el color desde los datos de la arista
        if color == "White":
            color = (1.0, 1.0, 1.0)
        elif color == "Black":
            color = (0.0, 0.0, 0.0)
        elif color == "Red":
            color = (1.0, 0.0, 0.0)
        elif color == "Green":
            color = (0.0, 1.0, 0.0)
        elif color == "Blue":
            color = (0.0, 0.0, 1.0)
        elif color == "Yellow":
            color = (1.0, 1.0, 0.0)
        elif color == "Magenta":
            color = (1.0, 0.0, 1.0)
        elif color == "Cyan":
            color = (0.0, 1.0, 1.0)
        elif color == "Gray":
            color = (0.5, 0.5, 0.5)
        elif color == "Orange":
            color = (1.0, 0.65, 0.0)
        mlab.plot3d(x, y, z, tube_radius=0.1, color=color)
    #Color nodos 
    if colorRNodo == "White":
        rgb_nr = (1.0, 1.0, 1.0)
    elif colorRNodo == "Black":
        rgb_nr = (0.0, 0.0, 0.0)
    elif colorRNodo == "Red":
        rgb_nr = (1.0, 0.0, 0.0)
    elif colorRNodo == "Green":
        rgb_nr = (0.0, 1.0, 0.0)
    elif colorRNodo == "Blue":
        rgb_nr = (0.0, 0.0, 1.0)
    elif colorRNodo == "Yellow":
        rgb_nr = (1.0, 1.0, 0.0)
    elif colorRNodo == "Magenta":
        rgb_nr = (1.0, 0.0, 1.0)
    elif colorRNodo == "Cyan":
        rgb_nr = (0.0, 1.0, 1.0)
    elif colorRNodo == "Gray":
        rgb_nr = (0.5, 0.5, 0.5)
    elif colorRNodo == "Orange":
        rgb_nr = (1.0, 0.65, 0.0)

    for node, (x, y, z) in posiciones.items():
        mlab.points3d(x, y, z, color=rgb_nr, scale_factor=1.0)
        etiqueta = f"{tr[node]}{num[node]}"
        mlab.text3d(x, y, z, etiqueta, scale=(0.4, 0.4, 0.4), color=(0.0, 0.0, 0.0))
    # Ajusta la vista 3D
    mlab.view(azimuth=0, elevation=90, distance=5)

    # Muestra la visualización interactiva
    mlab.show()

import tkinter as tk
from tkinter import ttk

Programa = tk.Tk()
Programa.title("BioNet")
Programa.iconbitmap("assets/BIONET.ico")

# Calcular el ancho y el alto de la ventana
window_width = Programa.winfo_screenwidth()  # Ancho de la mitad de la pantalla
window_height = Programa.winfo_screenheight()  # Altura de la pantalla

# Función para centrar la ventana en el lado izquierdo
def center_window(window, width, height):
    window.geometry(f"{width//2}x{height}")

center_window(Programa, window_width, window_height)

# Opciones de font: Fuente Arial, tamaño 12, negrita
label_font = ("sans-serif", 38, "bold")

# Opción de height: Altura del widget en caracteres de texto (3 líneas de alto)
label_height = 1

# Opción de width: Ancho del widget en caracteres de texto (20 caracteres de ancho)

Titulo = tk.Label(Programa, text="BioNet",font=label_font, height=label_height)
Titulo.place(x=200, y=10)

Logo = tk.PhotoImage(file='assets/BioNetLogo1.png')
LabelLogo=tk.Label(Programa,image=Logo)
LabelLogo.place(x=10, y=10)

#Seleccionar archivo
label_archivo = tk.Label(Programa, text="Selected file:")
label_archivo.place(x=200, y=100)

boton_seleccionar = tk.Button(Programa, text="Select file", command=BuscarArchivo)
boton_seleccionar.place(x=900, y=100)

label_id = tk.Label(Programa, text= "Enter the PDB file ID")
label_id.place(x=200, y=150)
 
id_pdb = tk.Entry(Programa)
id_pdb.place(x=400, y=150)
id_pdb.bind("<Return>", lambda event: IDPDB())

boton_buscar = tk.Button(Programa, text="Search file", command=DescargarArchivo)
boton_buscar.place(x=600, y=150)


# Crear un widget Text para mostrar el DataFrame
dataframe_text = tk.Text(Programa, height=25, width=130)
dataframe_text.place(x=10, y=400)

scroll_y = tk.Scrollbar(Programa, orient=tk.VERTICAL, command=dataframe_text.yview)
dataframe_text.configure(yscrollcommand=scroll_y.set)
scroll_y.pack(side=tk.RIGHT, fill=tk.Y)


###############################################################################################################
# Crear el widget Notebook
notebook = ttk.Notebook(Programa)
notebook.place(x=10, y=200)

################################################### Primera pestaña
tab1 = ttk.Frame(notebook, width=2000, height=150)
notebook.add(tab1, text="Plot Network")

label1 = tk.Label(tab1, text="Minimum threshold (Å):")
label1.place(x=10, y=40)

distancia = tk.Entry(tab1)
distancia.place(x=150, y=40)
distancia.bind("<Return>", lambda event: Distancia())

colores = ["White", "Black", "Red", "Green", "Blue", "Yellow", "Magenta", "Cyan", "Gray", "Orange"]
colores_e = ["White", "Black"]
LABELYN = ["Yes", "No"]

Color_Nodo = tk.Label(tab1, text="Select the node color:")
Color_Nodo.place(x=350, y=40)

Color_Arista = tk.Label(tab1, text="Select the edge color:")
Color_Arista.place(x=750, y=40)

opcion_Etiqueta = tk.Label(tab1, text="Do you want to add a label?")
opcion_Etiqueta.place(x=350, y=80)

Color_Etiqueta = tk.Label(tab1, text="Select the label color:")
Color_Etiqueta.place(x=750, y=80)

colores_n = ttk.Combobox(tab1, values=colores)
colores_n.place(x=550, y=40)

colores_a = ttk.Combobox(tab1, values=colores)
colores_a.place(x=970, y=40)

opcion_la = ttk.Combobox(tab1, values= LABELYN, state="readonly")
opcion_la.set("No")
opcion_la.place(x=550, y=80)

colores_l = ttk.Combobox(tab1, values=colores_e, state="disabled")
colores_l.set("White")
colores_l.place(x=970, y=80)

# Enlazar el evento <ComboboxSelected> al Combobox con la función seleccionar_elemento
colores_n.bind("<<ComboboxSelected>>", colornodo)
colores_a.bind("<<ComboboxSelected>>", colorarista)
colores_l.bind("<<ComboboxSelected>>", coloretiqueta)
opcion_la.bind("<<ComboboxSelected>>", lambda event: actualizar_estado_color())

#Boton Para constuir Red normal
btn_RED = tk.Button(tab1, text="Build Network", command=ConstruirRED)
btn_RED.place(x=150, y=80)

################################# Segunda pestaña
tab2 = ttk.Frame(notebook, width=2000, height=150)
notebook.add(tab2, text="Centrality properties")
#-------------------- Agregar la opcion para colorear con los tipos de residuos-----------------------------
propiedades_op = ["Degree", "Betweenness", "Closeness", "Clustering Coefficient", "Eigenvector", "Eigenvector Lamda2", "Eccentricity", "Stress"]

label1 = tk.Label(tab2, text="Minimum threshold (Å):")
label1.place(x=10, y=10)

distanciaP2 = tk.Entry(tab2)
distanciaP2.place(x=150, y=10)
distanciaP2.bind("<Return>", lambda event: Distancia2())

Label_Prop = tk.Label(tab2, text="Select a centrality property:")
Label_Prop.place(x=350, y=10)

propiedades_lista = ttk.Combobox(tab2, values=propiedades_op)
propiedades_lista.place(x=650, y=10)
propiedades_lista.bind("<<ComboboxSelected>>", propiedad_Funcion)

Label_Prop = tk.Label(tab2, text="Note: Display the centrality values before generating the graph.")
Label_Prop.place(x=1200, y=10)

boton_red_prop = tk.Button(tab2, text="Generate Network", command=Propiedad_Centralidad)
boton_red_prop.place(x=650, y=50)

boton_prop = tk.Button(tab2, text="Show selected centrality value", command=MostrarPropiedadesS)
boton_prop.place(x=850, y=5)

boton_centColDis = tk.Button(tab2, text="Show all centrality values", command=M_T_V)
boton_centColDis.place(x=850, y=50)

BotonGuardar = tk.Button(tab2, text= "Save data", command= Guardar)
BotonGuardar.place(x=850, y=100)
################################################# Tercer pestaña

tab3 = ttk.Frame(notebook, width=2000, height=150)
notebook.add(tab3, text="Distance Color")
####### Boton para generar la red de distancias predeterminada
boton_disR = tk.Button(tab3, text="Generate default distance network", command=DistanciasRed)
boton_disR.place(x=870, y=60)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Color y Rango de Aristas ~~~~~~~~~~~~~~~~~~~~~~~
coloresR = ["White", "Black", "Red", "Green", "Blue", "Yellow", "Magenta", "Cyan", "Gray", "Orange"]
rangos = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Entradas para rangos de valores ~~~~~~~~~~~~~~~~~~~~~~~~
distancia1_1 = ttk.Combobox(tab3, values=rangos, width=5)
distancia1_1.place(x=250, y=10)
distancia1_1.bind("<<ComboboxSelected>>", Rangos_Dis1_1)

distancia1_2 = ttk.Combobox(tab3, values=rangos, width=5)
distancia1_2.place(x=560, y=10)
distancia1_2.bind("<<ComboboxSelected>>", Rangos_Dis1_2)


distancia2_1 = ttk.Combobox(tab3, values=rangos, width=5)
distancia2_1.place(x=250, y=60)
distancia2_1.bind("<<ComboboxSelected>>", Rangos_Dis2_1)

distancia2_2 = ttk.Combobox(tab3, values=rangos, width=5)
distancia2_2.place(x=560, y=60)
distancia2_2.bind("<<ComboboxSelected>>", Rangos_Dis2_2)

distancia3_1 = ttk.Combobox(tab3, values=rangos, width=5)
distancia3_1.place(x=250, y=110)
distancia3_1.bind("<<ComboboxSelected>>", Rangos_Dis3_1)

distancia3_2 = ttk.Combobox(tab3, values=rangos, width=5)
distancia3_2.place(x=560, y=110)
distancia3_2.bind("<<ComboboxSelected>>", Rangos_Dis3_2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Etiquetas para rangos de valores ~~~~~~~~~~~~~~~~~~~~~~~~

Label_1min = tk.Label(tab3, text="Minimum threshold range 1 (Å):")
Label_1min.place(x=10, y=10)

Label_1max = tk.Label(tab3, text="Maximum threshold range 1 (Å):")
Label_1max.place(x=330, y=10)

Label_2min = tk.Label(tab3, text="Minimum threshold range 2 (Å):")
Label_2min.place(x=10, y=60)

Label_2max = tk.Label(tab3, text="Maximum threshold range 2 (Å):")
Label_2max.place(x=330, y=60)

Label_3min = tk.Label(tab3, text="Minimum threshold range 3 (Å):")
Label_3min.place(x=10, y=110)

Label_3max = tk.Label(tab3, text="Maximum threshold range 3 (Å):")
Label_3max.place(x=330, y=110)


#~~~~~~~~~~~~~~~~~~~~Etiquetas Colores aristas por rangos~~~~~~~~~~~~~~~~~~~~~~~~~~

Label_ColorRango1 = tk.Label(tab3, text="Color Range 1:")
Label_ColorRango1.place(x=640, y=10)

Label_ColorRango2 = tk.Label(tab3, text="Color Range 2:")
Label_ColorRango2.place(x=640, y=60)

Label_ColorRango3 = tk.Label(tab3, text="Color Range 3:")
Label_ColorRango3.place(x=640, y=110)

Label_ColorNodo = tk.Label(tab3, text="Node color:")
Label_ColorNodo.place(x=870, y=10)

#~~~~~~~~~~~~~~~~~~~~listas de Colores aristas por rangos~~~~~~~~~~~~~~~~~~~~~~~~~~
ColorRango1 = ttk.Combobox(tab3, values=coloresR, width=10)
ColorRango1.place(x=750, y=10)

ColorRango2 = ttk.Combobox(tab3, values=coloresR, width=10)
ColorRango2.place(x=750, y=60)

ColorRango3 = ttk.Combobox(tab3, values=coloresR, width=10)
ColorRango3.place(x=750, y=110)

ColorNodo = ttk.Combobox(tab3, values=coloresR)
ColorNodo.place(x=1030, y=10)

ColorRango1.bind("<<ComboboxSelected>>", Colores_Rango1)
ColorRango2.bind("<<ComboboxSelected>>", Colores_Rango2)
ColorRango3.bind("<<ComboboxSelected>>", Colores_Rango3)
ColorNodo.bind("<<ComboboxSelected>>", Colores_Rango_Nodo)

#~~~~~~~~~~~~~~~~~~~~Boton red personalizada~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boton_disRpersonalizada = tk.Button(tab3, text="Generate custom distance network", command=DistanciasRedPersonalizada)
boton_disRpersonalizada.place(x=870, y=110)

#~~~~~~~~~~~~~~~~~Boton mostrar propiedades de centralidad
#boton_centColDis = tk.Button(tab3, text="Mostrar propiedades ", command=Prop_Color_Dis)
#boton_centColDis.place(x=9, y=200)


Programa.mainloop()

