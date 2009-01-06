import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.*;
import javax.swing.event.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.color.*;

import java.util.*;

/*
 * x Selecteer hoogte (dataset)
 * x Boolean voor shading
 * - smooth/flatshading
 * - diffuse color
 * - ambient color
 * - light position
 */ 

public class HeightplotSelectPanel extends ColormapSelectPanel  {
    private boolean somethingChangedInTheLightingModel = true;
    
    private JCheckBox shadingBox;
    private ColorButton diffuseButton;
    private ColorButton ambientButton;
    private ColorButton positionButton;
    private JRadioButton smoothShadingButton;
    private JRadioButton flatShadingButton;
            
    private JFrame frame;

    public HeightplotSelectPanel(int minColor, int maxColor, int colorCount, int colormap, JFrame frame) {
            super(minColor, maxColor, colorCount, colormap, frame);
            colormapPreviewPanel.setVisible(false);
            colormapSelectPanel.setVisible(false);            
            this.frame = frame;
            
            add(initShadingPanel());
    }

    public boolean isShadingEnabled() {
        return shadingBox.getModel().isSelected();
    }

    public void setShadingEnabled(boolean value) {
        shadingBox.getModel().setSelected(value);
    }
    
    public float[] getLightAmbientColor() {
        float[] c = ambientButton.getBackground().getColorComponents(null);
        return new float[] {c[0], c[1], c[2], 1.0f};
    }
    
    public float[] getLightDiffuseColor() {
        float[] c = diffuseButton.getBackground().getColorComponents(null);
        return new float[] {c[0], c[1], c[2], 1.0f};
    }
    
    public float[] getLightPosition() {
        float[] c = positionButton.getBackground().getColorComponents(null);
        return new float[] {c[0] - 0.5f, c[1] - 0.5f, c[2] - 0.5f, 0.0f};
    }
    
    public boolean getSomethingChangedInTheLightingModel() {
        return somethingChangedInTheLightingModel;
    }
    
    public void setSomethingChangedInTheLightingModel(boolean value) {
        somethingChangedInTheLightingModel = value;
    }
    
    public boolean getFlatShading() {
        return flatShadingButton.getModel().isSelected();
    }
    
    private JPanel initShadingPanel() {
        shadingBox = new JCheckBox("Enable shading", true);
        shadingBox.setActionCommand("SHADING");
        shadingBox.addChangeListener(new HeightplotChangeHandler(this));
        
        diffuseButton = new ColorButton(this, frame);
        diffuseButton.setBackground(Color.WHITE);
        
        ambientButton = new ColorButton(this, frame);
        ambientButton.setBackground(Color.WHITE);      
        
        positionButton = new ColorButton(this, frame);
        positionButton.setBackground(Color.GREEN);      
        
        flatShadingButton = new JRadioButton("Flat shading", false);
        flatShadingButton.addChangeListener(new HeightplotChangeHandler(this));
        
        smoothShadingButton = new JRadioButton("Smooth shading", true);
        smoothShadingButton.addChangeListener(new HeightplotChangeHandler(this));
        
        ButtonGroup bg = new ButtonGroup();
        bg.add(flatShadingButton);
        bg.add(smoothShadingButton);
        
        JPanel shadingPanel = new JPanel();
        shadingPanel.setBorder(new TitledBorder("Heightplot Options"));
        shadingPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        shadingPanel.setLayout(new SpringLayout());
	
	shadingPanel.add(shadingBox);
        shadingPanel.add(new JLabel("Diffuse color: "));
        shadingPanel.add(diffuseButton);
        shadingPanel.add(new JLabel("Ambient color: "));
        shadingPanel.add(ambientButton);
        shadingPanel.add(new JLabel("Light position: "));
        shadingPanel.add(positionButton);
        shadingPanel.add(flatShadingButton);
        shadingPanel.add(smoothShadingButton);
        
        SpringUtilities.makeCompactGrid(shadingPanel, 9, 1,  // rows, cols
                                                      6, 6,  // initX, initY
                                                      6, 6); // xPad, yPad
        return shadingPanel;
    }
    
    class HeightplotChangeHandler implements ChangeListener, ActionListener {
        private HeightplotSelectPanel caller;
        
        public HeightplotChangeHandler(HeightplotSelectPanel caller) {
            this.caller = caller;
        }
        
        public void stateChanged(ChangeEvent e) {
            setSomethingChangedInTheLightingModel(true);
        }
        
        public void actionPerformed(ActionEvent e) {
            setSomethingChangedInTheLightingModel(true);           
        }
    }
    
    class ColorButton extends JButton implements ActionListener,ChangeListener {
        private HeightplotSelectPanel caller;
        private ColorSelector colorselector;
        
        public ColorButton(HeightplotSelectPanel caller, JFrame frame) {
            this.caller = caller;
            addActionListener(this);
            colorselector = new ColorSelector(frame);
            colorselector.addChangeListener(this);
        }
        
        public void actionPerformed(ActionEvent e) {
            colorselector.setVisible(true);
        }
        
        public void stateChanged(ChangeEvent e) {
            if(e.getSource().getClass().getName().equals("java.awt.Color")) {
                setBackground((Color)e.getSource());
                caller.setSomethingChangedInTheLightingModel(true);
            }
        }
        
	public void paintComponent(Graphics gr) {
            float width = getWidth();
            float height = getHeight();

            Graphics2D g = (Graphics2D)gr;
            g.setColor(getBackground());
            g.fillRect(0, 0, (int)width, (int)height);
	}
    }
}
