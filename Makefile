rootPath = ./
include ./include.mk

all : ${binPath}/ortheus_core

#These are parameters for composing/compiling the graph

#Use this line to point at a graphviz tool
makeGraph = dot

#This parameter can be adjusted (>=1) to control the silent-state number vs. the
#number of transitions.
collapseCoefficient = 1000 

#The location of the model files used
modelPath = ./models
model = ${modelPath}/affineModel.sxpr
paramModel = ${modelPath}/affineModel.param
 
clean :
	rm -f ${binPath}/ortheus_core ${modelPath}/model.otra xyzModelC.c xyzModelC.h
	cd submodules/sonLib && ${MAKE} cP.clean

${binPath}/ortheus_core : *.c *.h xyzModelC.c xyzModelC.h ${basicLibsDependencies}
	mkdir -p ${binPath}
	${cxx} ${cflags} -I ${libPath} -o ${binPath}/ortheus_core *.c ${basicLibs} -lm

# stdin from </dev/null works around stray stdin read on OS/X that hangs backgroud
# jobs
xyzModelC.c xyzModelC.h : ${model} ${paramModel} TransducerComposer.py TransducerCompiler.py
	rm -f ${modelPath}/model.otra xyzModelC.c xyzModelC.h ${modelPath}/model.dot ${modelPath}/ortheusmodel.pdf
	${PYTHON} TransducerComposer.py ${model} ${modelPath}/temp.otra ${modelPath}/model.dot ${collapseCoefficient} </dev/null
	cat  ${paramModel} ${modelPath}/temp.otra > ${modelPath}/model.otra
	rm ${modelPath}/temp.otra
	${PYTHON} TransducerCompiler.py ${modelPath}/model.otra xyzModelC.c xyzModelC.h </dev/null
	#Use this line if you want the pretty picture
	#${makeGraph} ${modelPath}/model.dot -Tpdf > ${modelPath}/model.pdf

${basicLibsDependencies} :
	cd submodules/sonLib && ${MAKE} cP
  
test :
	#Running python allTests.py
	${PYTHON} allTests.py --testLength=SHORT --logDebug
