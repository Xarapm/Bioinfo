#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 01 OCTUBRE 23:59
## Se requiere la entrega de este script completado con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)

data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #para revisar las dimensiones de los datos 
head(data) #para ver las primeras lineas
tail(data) #para ver las últimas líneas 

# Hacemos un primer histograma para explorar los datos
hist(data) #para crear un histograma

# Transformamos los datos con un logaritmo 
data_log=log2(data) #para transformar los datos 
hist(data_log) #para hacer un nuevo histograma con los datos transformados
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
#cambia la distribución de la gráfica, a una más normal
#puede servir para normalizar los datos

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
boxplot(data_log) #para generar un boxplot de los datos transformados
boxplot(data_log,col=c("blue","blue","blue","orange","orange","orange")) #para cambiar el color de los WT y los K.O, los azules son los W.T, los naranjas son los K.O, esto permite diferenciarlos
boxplot(data_log,col=c("blue","blue","blue","orange","orange","orange"),main="GSE5583-boxplots",las=2) #añadimos main="name" para cambiar el nombres del gráfico y las=2 para cambiar los ejes del gráfico

# ¿Qué es un boxplot?
#un boxplot es un diagrama de cajas que representa una serie de números por sus cuartíles, ayuda a localizar la media, la mediana 

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
hc=hclust(as.dist(1-cor(data_log))) #para agrupar los datos en clusters
plot(hc, main="clustering") #para cambiar el nombre del gráfico
#la separación es correcta ya que el cluster de la muestra esta separado por los W.T y los K.O

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt<-data[,1:3]
ko<-data[,4:6]
class(wt) #para generar la matriz, class nos dice el tipo de datos que tenemos
head(wt) #para leer las primeras lineas de la matriz 
#para seleccionar únicamente las columnas con los datos correspondientes a los wt y los ko 
#los datos obtenidos han sido una matriz de las columnas de cada grupo wt y ko 

# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean=apply(wt,1,mean) #calcular el parametro que queramos sobre los datos que queramos, generando una variable nueva 
wt.mean #para ver la media correspondiente a cada gen
head(wt.mean) #para ver las primeras lineas

#usamos los mismos comandos que antes, para calcular la media correspondiente a cada gen, en el caso de los ko
ko.mean=apply(ko,1,mean)
ko.mean
head(ko.mean)

# ¿Cuál es la media más alta?
max(wt.mean)
max(ko.mean)
#estos comandos sirven para averiguar la media más alta en el caso de los wt y de los ko
#la media más alta en este caso es la de los ko 

# Ahora hacemos un scatter plot (gráfico de dispersión)
#scatter plor; gráfico que muestra una diagonal en comparaciçon con las medias de los dos grupos 
plot(ko.mean~wt.mean)

#comando para realizar el scatter plot
plot(ko.mean~wt.mean,xlab="wt",ylab="ko",main="GSE5583-Scatter")

# Añadir una línea diagonal con abline
abline(0,1,col="red") #para crear una liena que indiuqe la regresion del gráfico
#el 0 y 1 indican la dirección de la pendeidnte 
#el col red indica el color lde la linea 

# ¿Eres capaz de añadirle un grid?
grid()
#para añadir una cuadrícula al plot 

# Calculamos la diferencia entre las medias de las condiciones
diff.mean=wt.mean-ko.mean
#para crear una tabla con las diferencias de meadias 

# Hacemos un histograma de las diferencias de medias
#para crear un histograma de la diferencia de medias 
hist (diff.mean)

# Calculamos la significancia estadística con un t-test.
#esto sirve para averiguar si los datos corespondientes a cada gen son estadísticamente significativos 

# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
#porque sin transformar son datos más fiables 
# ¿Cuántas valores tiene cada considicón?
#cada condición (wt y ko), tiene 3 valores 
pvalue=NULL
tstat=NULL
for(i in 1 : nrow(data)) {#para cada gen
  x=wt[i,] #gene wt número i
  y=ko[i,] #gene ko número i
  #Hacemos el test 
  t=t.test(x, y)
  
  #Añadimos el p-value a la lista 
  pvalue[i]=t$p.value
  #Añadimos las estadisticas a la lista
  tstat[i]=t$statistic
}

 head(pvalue) #para comprobar que el bucle se ha hecho correctamente
 length(pvalue) #para saber el tamaño y que hemos hecho todos los calculos 
 #para cada gen hemos conseguido un pvalue; el cual nos dice si la diferencia de las medias es un valor significativo
# Ahora comprobamos que hemos hecho TODOS los cálculos

# Hacemos un histograma de los p-values.
 hist(pvalue) #para hacer el histograma de los pvalue
# ¿Qué pasa si le ponemos con una transformación de -log10?
hist(-log10(pvalue), col="blue")
#reducimos el gráfico a los valores que nos interesan, y lo cambaimos de color 

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean,-log10(pvalue),main="GSE5583-VOLCANO")
#volcano plot de la diferencia de medias en contra del logaritmo (-log10)

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff=2 
pvalue_cutoff=0.01
#para establecer los mín y máx 
abline(v=diff.mean_cutoff,col="blue",lwd=3)
#abline(v=-diff.mean_cutoff,col="red",lwd=3)
abline(h=-log10(pvalue_cutoff),col="green",lwd=3)
#para establecer los valores entre los rangos que hemos establecido antes 
#los genes significativos se van a encontrar en las dos cuadrículas superiores 

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean=abs(diff.mean)>= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])
#para guardar todos los genes cuyo valor sea = o mayor que las diferencias de medias 

# Ahora el filtro de p-value
filter_by_pvalue=pvalue<=pvalue_cutoff
dim(data[filter_by_pvalue, ])
#para exptraer todos los pvalues que sean menor o igual que el pvalue 
#dim; dimensiones, para buscar cuantos genes sobrepasan el pvalue 

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined=filter_by_diff.mean&filter_by_pvalue
filtered=data[filter_combined,]
dim(filtered)
head(filtered)
#combinarlos para quedarnos con los genes que compartan ambos valores (que pertenezcan al filtro combinado) y eliminar los que no tienen valores similares 

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean,-log10(pvalue), main="GSE5583-VOLCANO #2")
points(diff.mean[filter_combined],-log10(pvalue[filter_combined]),col="red")
#sobrexpresados van a estar en la parte negativa de la gráfica y los reprimidos vana estar en la parte positiva, ya que la diferencia ha sido los wt - ko y estos últimos tienen valores más altos

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?


# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
heatmap(filtered)
#diferencia entre sobreexpresados (rojos) y reprimidos, con una expresión más baja  (amarillos)

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap


# Lo guardamos en un archivo PDF


# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table(filtered,"GSE5583_DE.txt",sep="\t",
            quote=FALSE)
