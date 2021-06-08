#include "headers/mainwindow.h"
#include "ui_mainwindow.h"
#include "headers/solver.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    /* Тест */
    Solver solver(10, 10, 1., 1., 5, 5e-4, 101325., 2, 3, 4,
                  12., 4, 5, 6);
    solver.processTMA(0);
}

MainWindow::~MainWindow()
{
    delete ui;
}

