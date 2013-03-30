# Алгоритм обработки данных:
# 1) создаём в рабочей дирректории набор папок, каждая папка должна содержать данные одного съема,
# а именно дирректорию с сегментарными матрицами МЭГИ (имя дирректории храниться в переменной matrixDirectory,
# по умолчанию "mt") и файл ("SI.txt") со строкой, с которой будет проводиться корреляция ячеек матриц
# (изначально делалось для сравнения со стресс индексом Варикарда)
# !!! Внимание: количество файлов с матрицами одного съёма и количество чисел в файле "SI.txt" должно совпадать.
# 2) вызываем функцию calculateCorrelations(), в качестве аргумента передавая строку, содержающую путь к
# рабочей дирректории (по умолчанию используется текущая). В дирректории каждого съёма будет создана папка
# (её имя храниться в переменной grDirectory, по умолчанию "graph"), содержащая корреляционные матрицы по каждому съёму
# 3) вызываем функцию correlationEvaluation(), в качестве аргумента передавая строку, содержающую путь к
# рабочей дирректории (по умолчанию используется текущая). В рабочей дирректории будет создан набор файлов, содержащих
# медианы и дисперсии корреляций всех съёмов.

#======================================================================================
#  Переменные настройки скрипта
#======================================================================================

grDirectory <- "graph"      # Дирректория вывода
matrixDirectory <- "mt"     # Дирректория с матрицами МЭГИ
varicardData <- "SI.txt"    # Файл, содержащий данные, с которыми будет прводиться корреляция вариации ячеек матриц

#======================================================================================
#  Константы: имена столбцов, палитры
#======================================================================================

segs <- c('C1', 'C2-3', 'C4-5', 'C6-7', 'C8',
          'Th1', 'Th2','Th3','Th5','Th6','Th7','Th8','Th10','Th11','Th12',
          'L1', 'L2','L3','L4','L5',
          'S1', 'S2', 'S3', 'K')

funcs <- c('F1-1', 'F1-2', 'F1-3', 'F1-4', 'F1-5',
           'F2-1', 'F2-2', 'F2-3', 'F2-4', 'F2-5',
           'F3-1', 'F3-2', 'F3-3', 'F3-4', 'F3-5',
           'F4-1', 'F4-2', 'F4-3', 'F4-4', 'F4-5',
           'F5-1', 'F5-2', 'F5-3', 'F5-4', 'F5-5',
           'F6-1', 'F6-2', 'F6-3', 'F6-4', 'F6-5',
           'F7-1', 'F7-2', 'F7-3', 'F7-4', 'F7-5')

# Сегментарные комплексы
actTable <- list(c('C1', 'Th2', 'Th3', 'Th10', 'L4'),
                 c('C2-3', 'Th3', 'Th4', 'Th10', 'Th12', 'L4', 'L5'),
                 c('C4-5', 'Th4', 'Th5', 'Th12', 'L1', 'L5', 'S1'),
                 c('C6-7', 'Th6', 'L1', 'S1', 'S2'),
                 c('C8', 'Th6', 'Th7', 'L1', 'L2', 'S2', 'S3'),
                 c('Th1', 'Th8', 'L3', 'S4', 'S5'),
                 c('Th2','Th9', 'L3', 'L4', 'K'),
                 c('Th3', 'C1', 'C2-3', 'C4-5', 'Th10', 'Th11', 'Th12', 'L4', 'L5'),
                 c('Th5', 'C4-5', 'Th12', 'L5', 'S1'),
                 c('Th6', 'C6-7', 'L1', 'S1', 'S2'),
                 c('Th7', 'C6-7', 'C8', 'L2', 'S3'),
                 c('Th8','Th1', 'Th2', 'L3', 'L4', 'K'),
                 c('Th10', 'C1', 'Th1', 'Th2', 'Th3', 'L4'),
                 c('Th11', 'C2-3', 'Th3', 'Th4', 'L5'),
                 c('Th12', 'C4-5', 'Th5', 'L5', 'S1'),
                 c('L1', 'C6-7', 'Th6', 'S1', 'S2'),
                 c('L2', 'C6-7', 'C8', 'Th6', 'Th7', 'S2', 'S3'),
                 c('L3', 'C8', 'Th1', 'Th7', 'Th8', 'S3', 'K'),
                 c('L4', 'C1', 'Th2', 'Th9', 'Th10', 'K'),
                 c('L5','C2-3', 'Th3', 'Th4', 'Th11'),
                 c('S1', 'C4-5', 'Th5', 'Th12', 'L1'),
                 c('S2', 'C6-7', 'Th6', 'Th7', 'L1'),
                 c('S3', 'C6-7', 'C8', 'Th6', 'Th7', 'Th8', 'L1', 'L2', 'L3'),
                 c('K', 'Th1', 'Th2', 'Th8', 'Th9', 'L3', 'L4'))



# Палитра цветов сегментарной матрицы
matrPalette <- colorRampPalette(c("blue", "black", "red"), interpolate = "linear")
# Палитра цветов корреляционной матрицы
matrCorrPalette <- colorRampPalette(c("yellow", "green", "black", "blue", "red"), interpolate = "linear")

#======================================================================================
#  Функции вспомогательные
#======================================================================================

saveGraphMatrix <- function (filename, data, palette, name) {
  # Сохранение данных матрицы в png
  png(file = filename)
  heatmap(data, Rowv=NA, Colv=NA, col = palette, scale="column", main=name)
  dev.off()
}

setMatrixRCNames <- function (matrix) {
  # Присвоение имен строкам/столбцам матрицы
  rownames(matrix) <- funcs
  colnames(matrix) <- segs
  return(matrix)
}

reverseMatrix <- function(mat) {
  # Разворот значений столбцов чтобы они начинались с F1-1, а не F7-5
  for (col in 1:ncol(mat)) {
    mat[,col] <- rev(mat[,col])
  }
  return(mat)
}

readSegMatrix <- function(filesDir, lateral="left") {
  # Читает папку с файлами матриц и возвращает массив матриц
  # указанной латеральности
  if (lateral == "left") {
    str <- 1:35
  }
  else if (lateral == "right") {
    str <- 36:70
  }
  else {
    print ("Error: bad lateral")
    return (0)
  }
  files <- dir(filesDir)
  matrArray <- array(data = NA, dim = c(35, 24,length(files)))
  count <- 1 # Счётчик матриц
  for (file in files) {
    matFile <- paste(filesDir, file, sep = "/")
    # Чтение файла с матрицей
    common.table <- read.table(matFile, sep = ';')
    # Выделяем из общего файла нужную матрицу - левую или правую
    m <- data.matrix(common.table[str,1:24])
    # Разворот матрицы
    m <- reverseMatrix(m)
    # Заносим полученную матрицу в общий массив
    matrArray[,,count] <- m
    # Увеличиваем счётчик
    count <- count + 1
  }
  return (matrArray)
}

#======================================================================================
#  Функции анализа
#======================================================================================

getCellLateralCorrelation <- function(dir, cellFunc, cellSeg) {
  # Корреляция между правым и левым полушарием конкретной ячейки
  DirName <- paste(dir, matrixDirectory, sep = "/")
  # Создаём массивы, содержащие сегментарные матриц
  matrArrayL <- readSegMatrix(DirName, "left")
  matrArrayR <- readSegMatrix(DirName, "right")
  L <- matrArrayL[match(cellFunc, funcs), match(cellSeg, segs), ]
  R <- matrArrayR[match(cellFunc, funcs), match(cellSeg, segs), ]
  c <- cor(L, R, method = "spearman")
  return (c)
}

getBoxplotDifference <- function(firstEx, cellFunc, cellSeg) {
  # Дисперсный анализ для одной ячейки
  firstDirName <- paste(firstEx, matrixDirectory, sep = "/")
  # Создаём массивы, содержащие сегментарные матрицы первого эксперимента
  fMatrArrayL <- readSegMatrix(firstDirName, "left")
  fMatrArrayR <- readSegMatrix(firstDirName, "right")
  # получаем вектора всех вариаций нужной ячейки
  fL <- fMatrArrayL[match(cellFunc, funcs), match(cellSeg, segs), ]
  fR <- fMatrArrayR[match(cellFunc, funcs), match(cellSeg, segs), ]
  return ( list(fL, fR) )
}

plotCellVarience <- function(cellFunc, cellSeg) {
  files <- dir(".")
  va <- list()
  vaNames <- vector()
  # l <- list(f1=x,f2=y)
  # c(l, f3=u)
  for (f in files) {
    c <- getCellLateralCorrelation(f, cellFunc, cellSeg)
    l <- getBoxplotDifference(f, cellFunc, cellSeg)
    
  }
}

getCellHarmonicVariations <- function(cellFunc, cellSeg, leg="") {
  # Влияние воздействия на ячейку, отражающееся в ячейках гармоник
  dir.create(grDirectory)
  files <- dir(matrixDirectory)
  # Хэш сегментарных комплексов
  segActHash <- new.env(hash=T)
  c <- 1
  for (s in segs) {
    segActHash[[s]] <- actTable[[c]]
    c <- c + 1
  }
  # Создаём массивы, содержащие сегментарные матриц
  matrArrayL <- readSegMatrix(matrixDirectory, "left")
  matrArrayR <- readSegMatrix(matrixDirectory, "right")
  # получаем список гармоник сегмента
  segHarmLst <- segActHash[[cellSeg]]
  # матрица для значений по сегментарным комплексам
  HarmMatrxL <- array(data = NA, dim = c(length(files), length(segHarmLst)))
  HarmMatrxR <- array(data = NA, dim = c(length(files), length(segHarmLst)))
  colnames(HarmMatrxL) <- segHarmLst
  colnames(HarmMatrxR) <- segHarmLst
  # Заполняем марицу значениями
  for (date in 1:length(files)) {
    for (s in segHarmLst) {
      HarmMatrxL[date, s] <- matrArrayL[match(cellFunc, funcs), match(s, segs), date]
      HarmMatrxR[date, s] <- matrArrayR[match(cellFunc, funcs), match(s, segs), date]
    }
  }
  NameL <- paste(grDirectory, "HarmMatrxL.txt", sep = "/")
  NameR <- paste(grDirectory, "HarmMatrxR.txt", sep = "/")
  write.table(HarmMatrxL, file=NameL)
  write.table(HarmMatrxR, file=NameR)
  for (s in segHarmLst) {
    grapgName <- paste(s, "png", sep = ".")
    grapgName <- paste(grDirectory, grapgName, sep = "/")
    png(file = grapgName)
    # Получаем разброс значений для столбца
    g_range <- range(0, HarmMatrxL[, match(s, segHarmLst)], HarmMatrxR[, match(s, segHarmLst)])
    plot(HarmMatrxR[, s], type="o", col="blue", ylim=g_range, axes=FALSE, ann=FALSE)
    lines(HarmMatrxL[, s], type="o", pch=22, lty=2, col="red")
    axis(2, las=1, at=4*0:g_range[2])
    box()
    l = paste(leg, s, sep = " ")
    title(main=l)
    title(xlab="Time")
    title(ylab="Amp")
    dev.off()
  }
}

getCellVariations <- function(cellFunc, cellSeg) {
  # Вариации отдельных ячеек сегментарной матрицы
  dir.create(grDirectory)
  files <- dir(matrixDirectory)
  # Создаём массивы, содержащие сегментарные матрицы
  matrArrayL <- readSegMatrix(matrixDirectory, "left")
  matrArrayR <- readSegMatrix(matrixDirectory, "right")
  l <- matrArrayL[which(funcs == cellFunc), which(segs == cellSeg), ]
  r <- matrArrayR[which(funcs == cellFunc), which(segs == cellSeg), ]
  grapgName <- paste(grDirectory, "cellVarG.png", sep = "/")
  png(file = grapgName)
  plot(l, type="o", col="green")
  lines(r, type="o", col="blue")
  dev.off()
  grapgName <- paste(grDirectory, "cellDiff.png", sep = "/")
  png(file = grapgName)
  dif <- l - r
  plot(dif, type="o", col="green")
  dev.off()
  grapgName <- paste(grDirectory, "cellVarB.png", sep = "/")
  png(file = grapgName)
  boxplot(l,r)
  dev.off()
}


getCorrTheta <- function(funcVect=21:27, ret=FALSE) {
  # Степень синхронизации с тетта ритмом
  dir.create(grDirectory)
  files <- dir(matrixDirectory)
  #======================================================================================
  #  Чтение матриц МЭГИ и строки Варикарда
  #======================================================================================
  # Создаём массивы, содержащие сегментарные матрицы
  matrArrayL <- readSegMatrix(matrixDirectory, "left")
  matrArrayR <- readSegMatrix(matrixDirectory, "right")
  # Чтение данных Варикарда
  SI <- scan(varicardData, dec = ",", n = length(files))
  #======================================================================================
  #  Рассчёт
  #======================================================================================
  # Создаём матрицу, содержащую корреляции кусков сегментов левого и правого полушарий
  corVectMatr <- array(data = NA, dim = c(24,length(files)))
  for (mat in 1:length(files)){
    for (s in 1:24) {
      corVectMatr[s, mat]<- cor(matrArrayL[funcVect,s,mat], matrArrayR[funcVect,s,mat], method = "spearman")
    }
  }
  # Запись матрицы промежуточных корреляций в файл
  CorVectMatrName <- paste(grDirectory, "CorMatVect.txt", sep = "/")
  rownames(corVectMatr) <- segs
  write.table(corVectMatr, file=CorVectMatrName)
  # Создаём вектор корреляций между ИН варикарда и изменением корреляции левого/правого по сегментам
  corVect <- array(data = NA, dim = c(24, 1))
  rownames(corVect) <- segs
  for (s in 1:24) {
    corVect[s,1] <- cor(corVectMatr[s,], SI, method = "spearman")
  }
  CorVectName <- paste(grDirectory, "CorVect.txt", sep = "/")
  write.table(corVect, file=CorVectName)
  if (ret == TRUE) {
    return (corVect[,1])
  }
}

getSegmentDistribution <- function(workdir=".", funcVect=11:17) {
  # Перебор всех папок в указанной дирректории и рассчет корреляций
  # между степенью синхронизации левого и правого полушарий и SI Варикарда.
  # Каждая папка должна содержать подпапку "mt" с набором матриц в текстовом
  # виде и файл "SI.txt" c SI Варикарда.
  # Вывод происходит в подпапку "graph" каждой папки
  setwd(workdir)
  dirList <- list.dirs(recursive = FALSE)
  m <- array(data=NA, dim=c(length(dirList),24))
  count <- 1
  for (d in dirList) {
    setwd(d)
    m[count,] <- getCorrTheta(funcVect, ret=TRUE)
    setwd("../")
    count <- count + 1
  }
  colnames(m) <- segs
  png(file = "segDist.png")
  boxplot(m)
  dev.off()
}

getCorrMatrix <- function(corBorder=0.5) {
  #======================================================================================
  #  Чтение матриц МЭГИ и строки Варикарда
  #======================================================================================
  dir.create(grDirectory)
  files <- dir(matrixDirectory)
  # Создаём массивы для хранения сегментарных матриц
  matrArrayL <- array(data = NA, dim = c(35, 24,length(files)))
  matrArrayR <- array(data = NA, dim = c(35, 24,length(files)))
  count <- 1 # Счётчик матриц
  for (file in files) {
    matFile <- paste(matrixDirectory, file, sep = "/")
    # Чтение файла с матрицей
    common.table <- read.table(matFile, sep = ';')
    # Разбиваем прочитанную таблицу на 2 матрицы - сначала левую, потом правую
    mL <- data.matrix(common.table[1:35,1:24])
    mR <- data.matrix(common.table[36:70,1:24])
    # Разворот матриц
    mL <- reverseMatrix(mL)
    mR <- reverseMatrix(mR)
    # Построение графиков текущей матрицы
    grapgName <- paste(file, ".png", sep = "")
    grapgName <- paste(grDirectory, grapgName, sep = "/")
    png(file = grapgName)
    par(mfrow = c(2,1))
    hist(mL)
    hist(mR)
    dev.off()
    # Заносим полученные матрицы в общие массивы
    matrArrayL[,,count] <- mL
    matrArrayR[,,count] <- mR
    # Увеличиваем счётчик
    count <- count + 1
  }
  # Чтение данных Варикарда
  SI <- scan(varicardData, dec = ",", n = length(files))
  #======================================================================================
  #  Расчёт корреляционной матрицы
  #======================================================================================
  # Создаём матрицы корреляций
  corrMatL <- array(data = NA, dim = c(35, 24))
  corrMatR <- array(data = NA, dim = c(35, 24))
  # Присвоение имён столбцам/строкам
  corrMatL <- setMatrixRCNames(corrMatL)
  corrMatR <- setMatrixRCNames(corrMatR)
  # Рассчёт корреляционных матриц
  count <- 1 # Счётчик сильных корреляций
  # Набор векторов для формирования конечного data.frame максимальных корреляций
  corV <- vector()
  hemisphereV <- vector()
  funcNameV <- vector()
  segNameV <- vector()
  for (f in 1:35){
    for (s in 1:24) {
      corrMatL[f,s] <- cor(matrArrayL[f,s,], SI, method = "spearman")
      # Если корреляция сильная - заполняем соответствующие вектора
      if ( (corrMatL[f,s] > corBorder) || (corrMatL[f,s] < -corBorder) ) {
        hemisphereV[count] <- "L"
        funcNameV[count] <- rownames(corrMatL)[f]
        segNameV[count] <- colnames(corrMatL)[s]
        corV[count] <- corrMatL[f,s]
        count <- count + 1
      }
      corrMatR[f,s] <- cor(matrArrayR[f,s,], SI, method = "spearman")
      # Если корреляция сильная - заполняем соответствующие вектора
      if ( (corrMatR[f,s] > corBorder) || (corrMatR[f,s] < -corBorder) ) {
        hemisphereV[count] <- "R"
        funcNameV[count] <- rownames(corrMatR)[f]
        segNameV[count] <- colnames(corrMatR)[s]
        corV[count] <- corrMatR[f,s]
        count <- count + 1
      }
    }
  }
  # Создаём результирующий data.frame сильных корреляций
  strongCorrDF <- data.frame(hemisphereV, segNameV, funcNameV, corV)
  # Построение графиков корреляционных матриц
  grapgName <- paste(grDirectory, "corrL.png", sep = "/")
  saveGraphMatrix(grapgName, corrMatL, matrCorrPalette(256), 'Left hemisphere correlation')
  grapgName <- paste(grDirectory, "corrR.png", sep = "/")
  saveGraphMatrix(grapgName, corrMatR, matrCorrPalette(256), 'Right hemisphere correlation')
  # Вывод коллеляционных матриц в файлы
  corNameL <- paste(grDirectory, "corrL.txt", sep = "/")
  corNameR <- paste(grDirectory, "corrR.txt", sep = "/")
  write.table(corrMatL, file=corNameL)
  write.table(corrMatR, file=corNameR)
  # Вывод таблицы сильных корреляций в файл
  strongCorName <- paste(grDirectory, "corrTable.txt", sep = "/")
  write.table(strongCorrDF, file=strongCorName)
}

calculateCorrelations <- function(workdir=".", corBorder=0.5) {
  # Перебор всех папок в указанной дирректории и рассчет корреляций
  # между ячейками матриц МЭГИ и SI Варикарда.
  # Каждая папка должна содержать подпапку "mt" с набором матриц в текстовом
  # виде и файл "SI.txt" c SI Варикарда.
  # Вывод происходит в подпапку "graph" каждой папки
  setwd(workdir)
  dirList <- list.dirs(recursive = FALSE)
  for (d in dirList) {
    setwd(d)
    getCorrMatrix(corBorder)
    setwd("../")
  }
}

correlationEvaluation <- function(workdir=".") {
  # Функция перебора набора папок, содержащих в дирректории вывода
  # данные построенных корреляций и построение по ним обобщения
  setwd(workdir)
  dirList <- list.dirs(recursive = FALSE)
  # Создаём массивы для хранения корреляционных матриц
  corrArrayL <- array(data = NA, dim = c(35, 24,length(dirList)))
  corrArrayR <- array(data = NA, dim = c(35, 24,length(dirList)))
  count <- 1
  # Заполняем массивы матрицами корреляций
  for (d in dirList) {
    corLName <- paste(d, grDirectory, "corrL.txt", sep="/")
    corRName <- paste(d, grDirectory, "corrR.txt", sep="/")
    # Чтение файла с матрицей
    corrArrayL[,,count] <- as.matrix(read.table(corLName, head=TRUE))
    corrArrayR[,,count] <- as.matrix(read.table(corRName, head=TRUE))
    count <- count + 1
  }
  # Создаём матрицы медиан корреляций
  corrMedianMatL <- array(data = NA, dim = c(35, 24))
  corrMedianMatR <- array(data = NA, dim = c(35, 24))
  # Создаём матрицы дисперсий корреляций
  corrDispMatL <- array(data = NA, dim = c(35, 24))
  corrDispMatR <- array(data = NA, dim = c(35, 24))
  for (f in 1:35) {
    for (s in 1:24) {
      corrMedianMatL[f,s] <- median(corrArrayL[f,s,])
      corrMedianMatR[f,s] <- median(corrArrayR[f,s,])
      corrDispMatL[f,s] <- var(corrArrayL[f,s,])
      corrDispMatR[f,s] <- var(corrArrayR[f,s,])
    }
  }
  # Присваиваем имена столбцам/строкам матриц
  corrMedianMatL <- setMatrixRCNames(corrMedianMatL)
  corrMedianMatR <- setMatrixRCNames(corrMedianMatR)
  corrDispMatL <- setMatrixRCNames(corrDispMatL)
  corrDispMatR <- setMatrixRCNames(corrDispMatR)
  # Сохранение результата в графической и текстовой форме
  saveGraphMatrix("corMedL.png", corrMedianMatL, matrCorrPalette(256), 'Left hemisphere correlation median')
  saveGraphMatrix("corMedR.png", corrMedianMatR, matrCorrPalette(256), 'Right hemisphere correlation median')
  saveGraphMatrix("corrDispL.png", corrDispMatL, matrPalette(256), 'Left hemisphere correlation variance')
  saveGraphMatrix("corrDispR.png", corrDispMatR, matrPalette(256), 'Right hemisphere correlation variance')
  write.table(corrMedianMatL, file="corMedL.txt")
  write.table(corrMedianMatR, file="corMedR.txt")
  write.table(corrDispMatL, file="corrDispL.txt")
  write.table(corrDispMatR, file="corrDispR.txt")
}