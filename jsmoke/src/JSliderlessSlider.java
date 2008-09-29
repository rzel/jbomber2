import java.awt.*;
import javax.swing.*;
import javax.swing.plaf.basic.*;
import java.util.*;
import java.awt.GradientPaint;

class JSliderlessSlider extends JSlider
{
    private float[][] colors;
    private int banding;
    private int count;

  public JSliderlessSlider(BoundedRangeModel brm, float[][] colors, int count)
{
    super();
    this.colors = colors;
    this.count  = count;
}

  public JSliderlessSlider(BoundedRangeModel brm)
{
    super(brm);
}

  public JSliderlessSlider(int orientation)
{
    super(orientation);
}

  public JSliderlessSlider(int min, int max)
{
    super(min, max);
}

  public JSliderlessSlider(int min, int max, int value)
{
    super(min, max, value);
}

  public JSliderlessSlider(int orientation, int min, int max, int value)
{
    super(orientation, min, max, value);
}

  public void setColors(float[][] colors, int count) {
    this.colors = colors;
    repaint();
}

  public void setCount(int count) {
      this.count = count;
      repaint();
}

  public void paintComponent(Graphics gr)
{
    double width = getWidth();
    double height = getHeight();

    Graphics2D g = (Graphics2D)gr;
    g.setColor(getBackground());
    g.fillRect(0,0,(int)width,(int)height);

    float k = 0;
    for(int i = 0; i < width; ++i) {
    	double pos = ((double)i) / width;
    	pos *= count + 1 ; pos = (int)pos ; pos /= count;
    	int c = (int)((colors.length - 1) * pos);
    	g.setColor(new Color(colors[c][0], colors[c][1], colors[c][2]));
    	g.drawLine(i, 0, i, (int)height / 2);
}

    ((BasicSliderUI)getUI()).paintLabels(g);
    ((BasicSliderUI)getUI()).paintTicks(g);
}
}