import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget
from PyQt5.QtGui import QPainter, QPen
from PyQt5.QtCore import Qt, QPoint

class DrawingApp(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Drawing Application")
        self.setGeometry(100, 100, 800, 600)

        self.canvas = CanvasWidget(self)
        self.setCentralWidget(self.canvas)
        
class CanvasWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.setGeometry(0, 0, 800, 600)

        self.drawing = False
        self.last_point = QPoint()
        self.current_point = QPoint()
        self.shapes = []

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.drawing = True
            self.last_point = event.pos()
            self.current_point = event.pos()

    def mouseMoveEvent(self, event):
        if self.drawing:
            self.current_point = event.pos()
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.drawing = False
            shape = (self.last_point, self.current_point)
            self.shapes.append(shape)
            self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        pen = QPen()
        pen.setWidth(2)
        pen.setColor(Qt.black)
        painter.setPen(pen)

        for shape in self.shapes:
            painter.drawLine(shape[0], shape[1])

        if self.drawing:
            painter.drawLine(self.last_point, self.current_point)

def main():
    app = QApplication(sys.argv)
    window = DrawingApp()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()