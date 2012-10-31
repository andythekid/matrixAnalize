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

segs <- c('C1', 'C2', 'C4', 'C6', 'C7',
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
  heatmap(data, Rowv=NA, Colv=NA, col = palette, scale="column", margins=c(5,5), main=name)
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

#======================================================================================
#  Функции анализа
#======================================================================================

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