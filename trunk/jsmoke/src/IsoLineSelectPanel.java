import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

public class IsoLineSelectPanel extends JPanel implements ActionListener, ChangeListener {
    public static final int DATASET_RHO   = 1;
    public static final int DATASET_F     = 2;
    public static final int DATASET_V     = 4;
    public static final int DATASET_F_DIV = 8;
    public static final int DATASET_V_DIV = 16;

    public static final int SCALE_CLAMP = 1;
    public static final int SCALE_SCALE = 2;

    private UnboundedDoubleSpinnerModel minClampSelectSpinnerModel = new UnboundedDoubleSpinnerModel(0.0);
    private UnboundedDoubleSpinnerModel maxClampSelectSpinnerModel = new UnboundedDoubleSpinnerModel(1.0);            

    private UnboundedDoubleSpinnerModel minIsoValueSpinnerModel = new UnboundedDoubleSpinnerModel(0.0);
    private UnboundedDoubleSpinnerModel maxIsoValueSpinnerModel = new UnboundedDoubleSpinnerModel(1.0);             
    
    private ColorSelector colorselector;    
    
    private float[] isoLineColor = {0.0f, 0.0f, 0.0f};
    private JLabel  colorLabel   = new JLabel("Current color");
    
    private int dataset   = DATASET_RHO;
    private int scalemode = SCALE_SCALE;
    
    private boolean updateIsoTexture = true;
    
    public IsoLineSelectPanel(JFrame frame) {
	setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));        
        add(initDatasetSelectPanel());
	add(initScalingSelectPanel());
        add(initIsoValuesSelectPanel());
        add(initIsoColorSelectPanel(frame));
        
        minIsoValueSpinnerModel.setValue(new Double(0.6));
        maxIsoValueSpinnerModel.setValue(new Double(0.6));
    }

    public int getDataset() {
        return dataset;
    }

    public int getScalemode() {
        return scalemode;
    }

    public double getMinIsoValue() {
        return ((Double)minIsoValueSpinnerModel.getValue()).doubleValue();
    }

    public double getMaxIsoValue() {
        return ((Double)maxIsoValueSpinnerModel.getValue()).doubleValue();
    }
    
    public float[] getIsoLineColor() {
        return isoLineColor;
    }
    
    public void setUpdateIsoTexture(boolean value) {
        updateIsoTexture = value;
    }
    
    public boolean getUpdateIsoTexture() {
        boolean ret = updateIsoTexture;
        return ret;
    }
    
    private JPanel initDatasetSelectPanel() {
        JRadioButton rhoButton = new JRadioButton("rho");
        rhoButton.setSelected(true);
        rhoButton.addActionListener(this);
        rhoButton.setActionCommand("SET_DATASET_RHO");

        JRadioButton fButton = new JRadioButton("|f|");
        fButton.addActionListener(this);
        fButton.setActionCommand("SET_DATASET_F");

        JRadioButton vButton = new JRadioButton("|v|");
        vButton.addActionListener(this);
        vButton.setActionCommand("SET_DATASET_V");

        JRadioButton dvButton = new JRadioButton("div(v)");
        dvButton.addActionListener(this);
        dvButton.setActionCommand("SET_DATASET_V_DIV");

        JRadioButton dfButton = new JRadioButton("div(f)");
        dfButton.addActionListener(this);
        dfButton.setActionCommand("SET_DATASET_F_DIV");                

        ButtonGroup datasetSelectGroup = new ButtonGroup();
        datasetSelectGroup.add(rhoButton);
        datasetSelectGroup.add(fButton);
        datasetSelectGroup.add(vButton);
        datasetSelectGroup.add(dfButton);
        datasetSelectGroup.add(dvButton);

        JPanel datasetSelectPanel = new JPanel();
        datasetSelectPanel.setLayout(new BoxLayout(datasetSelectPanel, BoxLayout.X_AXIS));
        datasetSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        datasetSelectPanel.setBorder(new TitledBorder("Dataset"));
        datasetSelectPanel.add(rhoButton);
        datasetSelectPanel.add(fButton);
        datasetSelectPanel.add(vButton);
        datasetSelectPanel.add(dfButton);
        datasetSelectPanel.add(dvButton);
        return datasetSelectPanel;
    }
    
    private JPanel initScalingSelectPanel() {
        JCheckBox clampButton = new JCheckBox("clamping");
        clampButton.setMnemonic(KeyEvent.VK_C);
        clampButton.setActionCommand("SET_SCALE_CLAMP");
        clampButton.addActionListener(this);

        JCheckBox scaleButton = new JCheckBox("scaling");
        scaleButton.setMnemonic(KeyEvent.VK_S);
        scaleButton.setActionCommand("SET_SCALE_SCALE");
        scaleButton.setSelected(true);
        scaleButton.addActionListener(this);

        JPanel clampSelectPanel = new JPanel();
        clampSelectPanel.setLayout(new GridLayout(3, 2));
        clampSelectPanel.add(clampButton);
        clampSelectPanel.add(new JLabel());
        clampSelectPanel.add(new JLabel("Low: ", SwingConstants.RIGHT));
        clampSelectPanel.add(new JSpinner(minClampSelectSpinnerModel));
        clampSelectPanel.add(new JLabel("High: ", SwingConstants.RIGHT));
        clampSelectPanel.add(new JSpinner(maxClampSelectSpinnerModel));

        JPanel scaleSelectPanel = new JPanel();
        scaleSelectPanel.setLayout(new BoxLayout(scaleSelectPanel, BoxLayout.Y_AXIS));
        scaleButton.setAlignmentX(Component.LEFT_ALIGNMENT);
        clampSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        scaleSelectPanel.setBorder(new TitledBorder("Scaling:"));
        scaleSelectPanel.add(scaleButton);
        scaleSelectPanel.add(clampSelectPanel);
        return scaleSelectPanel;        
    }
    
    public JPanel initIsoValuesSelectPanel() {
        JPanel isoValueSelectPanel = new JPanel();
        isoValueSelectPanel.setLayout(new GridLayout(2, 2));
        
        isoValueSelectPanel.add(new JLabel("Low: ", SwingConstants.RIGHT));
        isoValueSelectPanel.add(new JSpinner(minIsoValueSpinnerModel));
        isoValueSelectPanel.add(new JLabel("High: ", SwingConstants.RIGHT));
        isoValueSelectPanel.add(new JSpinner(maxIsoValueSpinnerModel));
        
        isoValueSelectPanel.setBorder(new TitledBorder("Iso values:"));
        isoValueSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        return isoValueSelectPanel;
    }
    
    public JPanel initIsoColorSelectPanel(JFrame frame) {
        colorselector = new ColorSelector(frame);
        colorselector.addChangeListener(this);
        
        JPanel isoColorSelectPanel = new JPanel();
        
        JButton button = new JButton("Change Iso line color");
        button.setActionCommand("CHANGE_COLOR");
        button.addActionListener(this);
        
        isoColorSelectPanel.add(button);
        isoColorSelectPanel.add(colorLabel);
        
        isoColorSelectPanel.setBorder(new TitledBorder("Iso values:"));
        isoColorSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);        
        return isoColorSelectPanel;
    }
      
    public void actionPerformed(ActionEvent e) {
        if (e.getActionCommand().equals("SET_DATASET_RHO")) {
                dataset = DATASET_RHO;
        } else if (e.getActionCommand().equals("SET_DATASET_F")) {
                dataset = DATASET_F;
        } else if (e.getActionCommand().equals("SET_DATASET_V_DIV")) {
                dataset = DATASET_V_DIV;
        } else if (e.getActionCommand().equals("SET_DATASET_F_DIV")) {
                dataset = DATASET_F_DIV;
        } else if (e.getActionCommand().equals("SET_DATASET_V")) {
                dataset = DATASET_V;
        } else if (e.getActionCommand().equals("SET_SCALE_CLAMP")) {
                scalemode ^= SCALE_CLAMP;
        } else if (e.getActionCommand().equals("SET_SCALE_SCALE")) {
                scalemode ^= SCALE_SCALE;
        } else if (e.getActionCommand().equals("CHANGE_COLOR")) {
            colorselector.setVisible(true);
        }
        
        updateIsoTexture = true;        
    }
    
    public void stateChanged(ChangeEvent e) {
        if(e.getSource().getClass().getName().equals("java.awt.Color")) {
            isoLineColor = ((Color)e.getSource()).getColorComponents(isoLineColor);
            colorLabel.setForeground((Color)e.getSource());
            colorLabel.repaint();
            updateIsoTexture = true;
        }
    }
}
